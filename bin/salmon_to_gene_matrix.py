#!/usr/bin/env python3
"""Aggregate per-sample Salmon quant.sf (transcript-level) into a gene-level
count/TPM matrix -- a small, dependency-free tximport equivalent.

Transcript-to-gene mapping is resolved in this order:
  1. --gene-map, a TSV of transcript_id<TAB>gene_id, if given explicitly.
  2. Derived from --gff3: for each mRNA/transcript feature, its Parent= is the
     gene. GFF3s vary (see inspect_annotation.py for why); many published
     annotations, including the Pinus densiflora v1.0 GFF3 this pipeline was
     built against, give every CDS/exon/gene the SAME ID with no Parent= at
     all -- there, "transcript ID" and "gene ID" are the same string by
     construction, and the mapping degenerates to identity.
  3. If neither is given, or a quant.sf ID has no entry in the derived map,
     identity is assumed (transcript ID treated as its own gene) and reported,
     not silently applied.
"""
from __future__ import annotations

import argparse
import os
import re
import sys


def parse_gff3_map(path: str) -> dict[str, str]:
    mapping: dict[str, str] = {}
    parent_style = 0
    identity_style = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            if feature not in ("mRNA", "transcript"):
                continue
            attrs = fields[8]
            tid = gid = None
            for kv in attrs.split(";"):
                kv = kv.strip()
                if kv.startswith("ID="):
                    tid = kv[3:]
                elif kv.startswith("Parent="):
                    gid = kv[7:].split(",")[0]
            if tid and gid:
                mapping[tid] = gid
                parent_style += 1

    if parent_style == 0:
        # No mRNA/transcript features with Parent=, or none at all (the
        # Pinus densiflora v1.0 case: 44,233 genes, only 4,637 mRNA features,
        # none carrying Parent=). Fall back to gene IDs directly -- CDS
        # sequences in the shipped *.CDS.fa are keyed by the gene ID itself.
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != "gene":
                    continue
                for kv in fields[8].split(";"):
                    kv = kv.strip()
                    if kv.startswith("ID="):
                        gid = kv[3:]
                        mapping[gid] = gid
                        identity_style += 1
    return mapping, parent_style, identity_style


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quant", nargs="+", required=True,
                    help="sample_id=path/to/quant.sf pairs")
    ap.add_argument("--gene-map", default=None)
    ap.add_argument("--gff3", default=None)
    ap.add_argument("--out-counts", required=True)
    ap.add_argument("--out-tpm", required=True)
    args = ap.parse_args()

    mapping: dict[str, str] = {}
    if args.gene_map:
        with open(args.gene_map) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    mapping[parts[0]] = parts[1]
        sys.stderr.write(f"INFO  ~ [salmon] gene map: {len(mapping):,} entries from {args.gene_map}\n")
    elif args.gff3:
        mapping, n_parent, n_identity = parse_gff3_map(args.gff3)
        if n_parent:
            sys.stderr.write(
                f"INFO  ~ [salmon] gene map: {n_parent:,} transcript->gene pairs "
                f"from Parent= in {args.gff3}\n"
            )
        if n_identity:
            sys.stderr.write(
                f"WARN  ~ [salmon] {args.gff3} has no usable mRNA Parent= links; "
                f"treating {n_identity:,} gene IDs as their own transcript "
                f"(matches the shipped CDS FASTA convention for this annotation)\n"
            )

    samples: list[str] = []
    counts: dict[str, dict[str, float]] = {}
    tpm: dict[str, dict[str, float]] = {}
    identity_fallback_hits = 0

    for pair in args.quant:
        sample, _, path = pair.partition("=")
        if not path:
            sys.stderr.write(f"ERROR ~ [salmon] malformed --quant entry {pair!r}, expected sample=path\n")
            return 1
        samples.append(sample)
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            idx = {name: i for i, name in enumerate(header)}
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                tid = fields[idx["Name"]]
                gid = mapping.get(tid)
                if gid is None:
                    gid = tid
                    identity_fallback_hits += 1
                counts.setdefault(gid, {})[sample] = counts.get(gid, {}).get(sample, 0.0) + float(fields[idx["NumReads"]])
                tpm.setdefault(gid, {})[sample] = tpm.get(gid, {}).get(sample, 0.0) + float(fields[idx["TPM"]])

    if identity_fallback_hits and mapping:
        sys.stderr.write(
            f"WARN  ~ [salmon] {identity_fallback_hits:,} transcript IDs across all "
            f"samples were not found in the gene map; treated as their own gene\n"
        )

    genes = sorted(counts)
    for path, table, label in [(args.out_counts, counts, "NumReads"),
                                (args.out_tpm, tpm, "TPM")]:
        with open(path, "w") as fh:
            fh.write("gene_id\t" + "\t".join(samples) + "\n")
            for gene in genes:
                row = table.get(gene, {})
                fh.write(gene + "\t" + "\t".join(f"{row.get(s, 0.0):.3f}" for s in samples) + "\n")
        sys.stderr.write(f"INFO  ~ [salmon] wrote {label} matrix: {len(genes):,} genes x {len(samples)} samples -> {path}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
