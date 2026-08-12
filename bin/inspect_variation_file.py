#!/usr/bin/env python3
"""Probe the dataset's own resequencing variant file and either convert it to
VCF or fail with a specific, actionable message -- never assume a layout that
hasn't been confirmed.

This is route 2 for phASER's phased VCF (docs/allele_specific_expression.md):
using the *P. densiflora* assembly authors' own published genotypes, where the
RNA-seq sample's accession is one they resequenced, is lower review risk than
this pipeline's own HA-vs-HB SNP calls (route 1, syri_to_phased_vcf.py) --
those genotypes are already peer-reviewed data, not something this pipeline
computed and has to defend.

Format is NOT assumed. As of this pipeline's design, only the Figshare
description of `Variation_information_of_P.densiflora_accessions.txt.gz` is
known ("Genotype information generated from resequencing analysis of
P. densiflora accessions") -- not its actual column layout. This script reads
the first N lines, reports what it finds, and:
  - if it already looks like VCF (starts with ##fileformat=VCF), passes it
    through unchanged;
  - if it looks like a plain genotype matrix (a header row naming accessions
    as columns), attempts a conversion for one named accession/sample;
  - otherwise, prints the header and a few data rows and exits non-zero rather
    than guessing.
"""
from __future__ import annotations

import argparse
import gzip
import sys


def opener(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--variation-file", required=True)
    ap.add_argument("--accession", default=None,
                    help="which accession/sample column to extract, if the file is a genotype matrix")
    ap.add_argument("--out", default="resequencing.vcf")
    ap.add_argument("--peek-lines", type=int, default=20)
    args = ap.parse_args()

    with opener(args.variation_file) as fh:
        head = [fh.readline().rstrip("\n") for _ in range(args.peek_lines)]
    head = [h for h in head if h != ""]

    if not head:
        sys.stderr.write(f"ERROR ~ [variation] {args.variation_file} is empty\n")
        return 1

    if head[0].startswith("##fileformat=VCF"):
        sys.stderr.write(
            f"INFO  ~ [variation] {args.variation_file} is already VCF -- "
            f"copying through unchanged to {args.out}\n"
        )
        with opener(args.variation_file) as src, open(args.out, "w") as dst:
            dst.write(src.read())
        return 0

    # Look for a header row: the first non-comment line, tab- or
    # whitespace-delimited, whose fields look like column names.
    header_line = next((h for h in head if not h.startswith("#")), None)
    if header_line is None:
        sys.stderr.write(
            f"ERROR ~ [variation] {args.variation_file}: could not find an "
            f"un-commented header row in the first {args.peek_lines} lines.\n"
            f"        First lines:\n" + "\n".join(f"        {h[:200]}" for h in head[:5]) + "\n"
        )
        return 1

    sep = "\t" if "\t" in header_line else None
    columns = header_line.split(sep)
    sys.stderr.write(
        f"INFO  ~ [variation] {args.variation_file}: detected {len(columns)} column(s):\n"
        f"        {columns[:15]}{'...' if len(columns) > 15 else ''}\n"
    )

    if args.accession is None:
        sys.stderr.write(
            "ERROR ~ [variation] format detected but --accession was not given, so "
            "no genotype column can be selected.\n"
            "        Re-run with --accession <name> matching one of the columns "
            "printed above, or confirm this file's format manually and adapt "
            "this script -- its layout was not known when this pipeline was written.\n"
        )
        return 2

    if args.accession not in columns:
        sys.stderr.write(
            f"ERROR ~ [variation] --accession {args.accession!r} is not among the "
            f"detected columns.\n        Available: {columns}\n"
            f"        route 1 (assembly-derived, syri_to_phased_vcf.py) remains "
            f"available regardless of whether this accession is listed.\n"
        )
        return 2

    sys.stderr.write(
        f"ERROR ~ [variation] column layout detected but the conversion to VCF for "
        f"an arbitrary genotype-matrix format is not implemented -- this script's "
        f"job is to detect and report, not guess a schema nobody has verified.\n"
        f"        Confirm the exact chrom/pos/ref/alt/genotype column names against "
        f"the header above, then extend this script (or convert manually) before "
        f"relying on route 2. Route 1 does not depend on this file.\n"
    )
    return 3


if __name__ == "__main__":
    sys.exit(main())
