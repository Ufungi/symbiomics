#!/usr/bin/env python3
"""Merge per-sample count files into one matrix, plus a QC summary.

Handles featureCounts and HTSeq output. Gene order is taken from the first file
and every subsequent file is required to carry the same gene set -- a silent
mismatch would produce a matrix where a column is offset against the rows,
which is the kind of error that survives all the way into a figure.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

# HTSeq writes its summary rows inline with the counts.
HTSEQ_SPECIAL = {
    "__no_feature", "__ambiguous", "__too_low_aQual",
    "__not_aligned", "__alignment_not_unique",
}


def read_featurecounts(path: str):
    counts, sample = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "Geneid":
                sample = os.path.basename(fields[-1])
                for suffix in (".sorted.bam", ".bam"):
                    if sample.endswith(suffix):
                        sample = sample[: -len(suffix)]
                continue
            counts[fields[0]] = int(float(fields[-1]))
    return sample, counts, {}


def read_htseq(path: str):
    counts, special = {}, {}
    sample = os.path.basename(path)
    for suffix in (".htseq.tsv", ".tsv", ".txt"):
        if sample.endswith(suffix):
            sample = sample[: -len(suffix)]
            break
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            name, value = fields[0], int(float(fields[1]))
            if name in HTSEQ_SPECIAL:
                special[name.lstrip("_")] = value
            else:
                counts[name] = value
    return sample, counts, special


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tool", required=True, choices=["featurecounts", "htseq"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--stats", required=True)
    ap.add_argument("--strandedness", default=None,
                    help="TSV: sample_id<TAB>strandedness[<TAB>conflict]")
    ap.add_argument("files", nargs="+")
    args = ap.parse_args()

    reader = read_featurecounts if args.tool == "featurecounts" else read_htseq

    order: list[str] = []
    per_sample: dict[str, dict[str, int]] = {}
    specials: dict[str, dict[str, int]] = {}

    for path in sorted(args.files):
        sample, counts, special = reader(path)
        if not sample:
            sample = os.path.basename(path).split(".")[0]
        if sample in per_sample:
            sys.stderr.write(f"ERROR ~ [counts] duplicate sample column {sample!r}\n")
            return 1
        if not order:
            order = list(counts.keys())
        elif set(counts.keys()) != set(order):
            missing = len(set(order) - set(counts.keys()))
            extra = len(set(counts.keys()) - set(order))
            sys.stderr.write(
                f"ERROR ~ [counts] {path}: gene set differs from the first file "
                f"({missing} missing, {extra} extra).\n"
                f"        All samples must be counted against the same annotation.\n"
            )
            return 1
        per_sample[sample] = counts
        specials[sample] = special

    samples = sorted(per_sample)
    with open(args.out, "w") as fh:
        fh.write("gene_id\t" + "\t".join(samples) + "\n")
        for gene in order:
            fh.write(gene + "\t" + "\t".join(str(per_sample[s][gene]) for s in samples) + "\n")

    # Strandedness comes in as a 2-column TSV covering EVERY sample, built from
    # the effective per-sample value. Reading the inference JSONs instead would
    # only describe the samples that were probed, leaving declared ones blank.
    strand_map = {}
    if args.strandedness and os.path.exists(args.strandedness):
        with open(args.strandedness) as fh:
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 2 and fields[0] != "sample":
                    strand_map[fields[0]] = {
                        "effective": fields[1],
                        "conflict": len(fields) > 2 and fields[2] == "yes",
                    }

    with open(args.stats, "w") as fh:
        header = ["sample", "assigned", "genes_detected", "strandedness", "strand_conflict"]
        if args.tool == "htseq":
            header += sorted({k for s in specials.values() for k in s})
        fh.write("\t".join(header) + "\n")
        for sample in samples:
            counts = per_sample[sample]
            assigned = sum(counts.values())
            detected = sum(1 for v in counts.values() if v > 0)
            info = strand_map.get(sample, {})
            row = [sample, str(assigned), str(detected),
                   info.get("effective", "NA"),
                   "yes" if info.get("conflict") else "no"]
            if args.tool == "htseq":
                row += [str(specials[sample].get(k, 0)) for k in header[5:]]
            fh.write("\t".join(row) + "\n")

    # A run that mixes forward and reverse libraries must not silently produce
    # one matrix -- that is a DE analysis that is wrong in a way nobody sees.
    calls = {v.get("effective") for v in strand_map.values()
             if v.get("effective") in ("forward", "reverse")}
    if len(calls) > 1:
        sys.stderr.write(
            f"WARN  ~ [counts] samples disagree on strandedness ({sorted(calls)}). "
            f"Mixing them in one differential-expression analysis is not valid.\n"
        )

    sys.stderr.write(
        f"INFO  ~ [counts] {len(order)} features x {len(samples)} samples -> {args.out}\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
