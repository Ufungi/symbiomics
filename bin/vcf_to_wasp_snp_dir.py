#!/usr/bin/env python3
"""Convert a phased VCF (syri_to_phased_vcf.py's output, or any biallelic SNP
VCF) into WASP's `--snp_dir` layout: one file per chromosome named
`<chrom>.snps.txt`, three whitespace-separated columns (position, ref, alt),
one row per SNP -- the simpler of WASP's two supported input formats (the
other needs a trio of HDF5 files that gains nothing for a single individual's
two known haplotypes).
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys


def opener(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    handles: dict[str, "TextIO"] = {}
    n = 0
    with opener(args.vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            chrom, pos, _id, ref, alt = fields[:5]
            if len(ref) != 1 or len(alt) != 1:
                continue  # WASP's text format is SNP-only
            if chrom not in handles:
                handles[chrom] = open(os.path.join(args.out_dir, f"{chrom}.snps.txt"), "w")
            handles[chrom].write(f"{pos}\t{ref}\t{alt}\n")
            n += 1

    for handle in handles.values():
        handle.close()

    sys.stderr.write(
        f"INFO  ~ [wasp] {n:,} SNPs across {len(handles)} chromosome file(s) -> {args.out_dir}\n"
    )
    if n == 0:
        sys.stderr.write(
            "WARN  ~ [wasp] 0 SNPs written -- WASP filtering will have nothing to check "
            "and every read will pass through unfiltered. Check the input VCF.\n"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
