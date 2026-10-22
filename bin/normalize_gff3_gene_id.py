#!/usr/bin/env python3
"""Ensure count_feature lines carry a gene-level ``gene_id`` attribute.

HTSeq and featureCounts group reads by whatever attribute sits on the feature
line being counted. The published model is ``--type exon --idattr gene_id``:
exons are the intervals, gene_id is the key, so every exon of one gene folds
into a single row. Neither tool walks a GFF3 ``Parent=`` hierarchy.

Many GFF3 files (funannotate, BRAKER3 via gffread, standard GFF3) put the gene
identity only on the gene line and link children through Parent:

    gene  ID=g1
    mRNA  ID=g1.t1;Parent=g1
    exon  ID=g1.t1.exon1;Parent=g1.t1          <- no gene_id

Counting that file with ``--idattr gene_id`` falls back to the exon's own
``ID`` and emits one row per exon. This script resolves each count_feature
line back to its gene and appends ``gene_id=...`` so the official flags work
without changing count_feature away from exon.

Already-normalized files (GTF, Ensembl GFF3 with gene_id on exons) pass
through. Lines that cannot be tied to a gene are an error when they are too
many to ignore -- a silent fallback is how the wrong row count ships.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import sys

# Coverage required before we refuse to emit an annotation. Mirrors
# inspect_annotation.py: a handful of odd lines must not fail the run, but
# counting the wrong feature level must.
MIN_COVERAGE = 0.99


def opener(path: str, mode: str = "rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_attrs(field: str) -> dict:
    out = {}
    # GTF: key "value";   GFF3: key=value;
    if "=" in field and not re.search(r'\w\s+"', field):
        for kv in field.split(";"):
            kv = kv.strip()
            if "=" in kv:
                key, _, value = kv.partition("=")
                out[key.strip()] = value.strip()
    else:
        for kv in field.split(";"):
            kv = kv.strip()
            match = re.match(r'^(\S+)\s+"?([^"]*)"?$', kv)
            if match:
                out[match.group(1)] = match.group(2)
    return out


def format_attrs(attrs: dict, style: str) -> str:
    if style == "gtf":
        return " ".join(f'{k} "{v}";' for k, v in attrs.items())
    return ";".join(f"{k}={v}" for k, v in attrs.items()) + ";"


def attr_style(attrs: dict) -> str:
    """GTF when the line uses transcript_id=... and has no GFF3 Parent chain."""
    if "transcript_id" in attrs and "Parent" not in attrs:
        return "gtf"
    return "gff3"


class FeatureIndex:
    """ID -> record, plus the gene IDs declared as type=gene."""

    def __init__(self) -> None:
        self.by_id: dict[str, dict] = {}
        self.gene_ids: set[str] = set()

    def add(self, fid: str, ftype: str, parent: str | None, gene_id: str | None) -> None:
        rec = self.by_id.get(fid)
        if rec is None:
            self.by_id[fid] = {
                "type": ftype,
                "parent": parent or "",
                "gene_id": gene_id or "",
            }
        else:
            if parent and not rec["parent"]:
                rec["parent"] = parent
            if gene_id and not rec["gene_id"]:
                rec["gene_id"] = gene_id
            # Prefer the gene typing when an ID is reused across levels.
            if ftype == "gene":
                rec["type"] = "gene"
        if ftype == "gene":
            self.gene_ids.add(fid)

    def resolve_gene(self, fid: str, seen: set[str] | None = None) -> str | None:
        """Walk Parent= up to the gene; also accept Pinus-style ID reuse."""
        if not fid:
            return None
        if seen is None:
            seen = set()
        if fid in seen:
            return None
        seen.add(fid)

        if fid in self.gene_ids:
            return fid

        rec = self.by_id.get(fid)
        if rec is None:
            return None
        if rec["gene_id"]:
            return rec["gene_id"]
        if rec["type"] == "gene":
            return fid
        if rec["parent"]:
            for parent in rec["parent"].split(","):
                parent = parent.strip()
                if not parent:
                    continue
                gene = self.resolve_gene(parent, seen)
                if gene:
                    return gene
        return None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--out", default="annotation.normalized.gff3")
    ap.add_argument("--want-feature", default="exon",
                    help="feature type that counters will read (params.count_feature)")
    ap.add_argument("--report", default=None,
                    help="optional JSON summary path")
    ap.add_argument("--max-lines", type=int, default=0,
                    help="0 = read the whole file")
    args = ap.parse_args()

    # --- pass 1: index every feature that carries ID=/gene_id= -------------
    index = FeatureIndex()
    raw_lines: list[str] = []
    # (fields, attrs) for data lines; None for headers/comments
    parsed: list[tuple[list[str], dict] | None] = []
    n_want = 0
    n_want_with_gene_id = 0

    with opener(args.annotation) as fh:
        for i, line in enumerate(fh):
            if args.max_lines and i >= args.max_lines:
                break
            if not line.strip() or line.startswith("#"):
                raw_lines.append(line)
                parsed.append(None)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                raw_lines.append(line)
                parsed.append(None)
                continue
            attrs = parse_attrs(fields[8])
            raw_lines.append(line)
            parsed.append((fields, attrs))

            fid = attrs.get("ID", "")
            if fid:
                index.add(fid, fields[2], attrs.get("Parent"), attrs.get("gene_id"))
            elif attrs.get("gene_id"):
                # GTF exons often lack ID; gene_id itself is the key.
                index.add(attrs["gene_id"], fields[2],
                          attrs.get("Parent"), attrs["gene_id"])

            if fields[2] == args.want_feature:
                n_want += 1
                if attrs.get("gene_id"):
                    n_want_with_gene_id += 1

    if not raw_lines:
        sys.stderr.write(
            f"ERROR ~ [normalize] {args.annotation} has no readable lines.\n"
        )
        return 1

    if n_want == 0:
        # Nothing to normalize; let INSPECT_ANNOTATION report the real issue.
        sys.stderr.write(
            f"WARN  ~ [normalize] no '{args.want_feature}' features in "
            f"{args.annotation}; passing the file through unchanged\n"
        )
        with opener(args.out, "wt") as out:
            out.writelines(raw_lines)
        _write_report(args, n_want, n_want_with_gene_id, 0, 0, [], passthrough=True)
        return 0

    # --- pass 2: append gene_id on count_feature lines that lack it --------
    out_lines: list[str] = []
    n_added = 0
    n_failed = 0
    failed_ids: list[str] = []

    for entry, line in zip(parsed, raw_lines):
        if entry is None:
            out_lines.append(line)
            continue
        fields, attrs = entry
        if fields[2] != args.want_feature or attrs.get("gene_id"):
            out_lines.append(line)
            continue

        gene = None
        # Pinus densiflora v1.0: exon lines repeat the gene's own ID.
        fid = attrs.get("ID", "")
        if fid and fid in index.gene_ids:
            gene = fid
        if gene is None and fid:
            gene = index.resolve_gene(fid)
        if gene is None:
            for parent in (attrs.get("Parent") or "").split(","):
                parent = parent.strip()
                if parent:
                    gene = index.resolve_gene(parent)
                    if gene:
                        break

        if gene is None:
            n_failed += 1
            if len(failed_ids) < 10:
                failed_ids.append(fid or attrs.get("Parent") or fields[0])
            out_lines.append(line)
            continue

        attrs["gene_id"] = gene
        fields[8] = format_attrs(attrs, attr_style(attrs))
        out_lines.append("\t".join(fields) + "\n")
        n_added += 1

    coverage = (n_want_with_gene_id + n_added) / n_want if n_want else 0.0
    if coverage < MIN_COVERAGE:
        sys.stderr.write(
            f"ERROR ~ [normalize] only {coverage:.1%} of '{args.want_feature}' "
            f"lines in {args.annotation} carry gene_id after Parent resolution\n"
            f"        added={n_added} already={n_want_with_gene_id} "
            f"failed={n_failed} total={n_want}\n"
            f"        first unresolved: {', '.join(failed_ids) or '(none)'}\n"
            f"        Counting would group at the wrong level. Fix the GFF3 "
            f"hierarchy or set count_attribute explicitly.\n"
        )
        return 1

    with opener(args.out, "wt") as out:
        out.writelines(out_lines)

    sys.stderr.write(
        f"INFO  ~ [normalize] '{args.want_feature}' lines with gene_id: "
        f"{n_want_with_gene_id + n_added:,}/{n_want:,} "
        f"(added {n_added:,}, unchanged {n_want_with_gene_id:,})\n"
    )
    if n_failed:
        sys.stderr.write(
            f"WARN  ~ [normalize] {n_failed:,} line(s) lacked gene_id and "
            f"could not be resolved; coverage still {coverage:.1%}\n"
        )

    _write_report(args, n_want, n_want_with_gene_id + n_added, n_added,
                  n_failed, failed_ids, passthrough=False)
    return 0


def _write_report(args, n_want, n_covered, n_added, n_failed, failed_ids,
                  passthrough: bool) -> None:
    if not args.report:
        return
    payload = {
        "annotation": args.annotation,
        "want_feature": args.want_feature,
        "n_want_feature": n_want,
        "n_with_gene_id": n_covered,
        "n_added": n_added,
        "n_failed": n_failed,
        "unresolved_examples": failed_ids,
        "passthrough": passthrough,
    }
    with open(args.report, "w") as fh:
        json.dump(payload, fh, indent=2)


if __name__ == "__main__":
    sys.exit(main())
