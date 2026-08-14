#!/usr/bin/env python3
"""Turn SyRI's syri.out SNP rows into a phased VCF for phASER.

syri.out is a 12-column TSV (SyRI docs, confirmed format):
  1 ref_chr  2 ref_start  3 ref_end  4 ref_seq  5 query_seq  6 query_chr
  7 query_start  8 query_end  9 ID  10 parent_ID  11 annotation_type
  12 copy_status

For SNP rows, columns 4/5 are the literal REF (HA) and ALT (HB) alleles at
ref_start on ref_chr. Phase is known by construction, not statistically
inferred: HA and HB are already the two separated haplotypes of this
individual, so every SNP is unambiguously HA=0, HB=1 -- GT is always 0|1, never
0/1 (unphased) or 1|0. This is the property that makes the assembly-alignment
route cheaper and more defensible than aligning RNA-seq to one diploid
reference and statistically phasing heterozygous calls.

Restricting to gene-overlapping positions (--genes) is optional but recommended
-- phASER only needs sites inside transcribed regions, and it shrinks the VCF a
lot on a 21.7 Gb genome.
"""
from __future__ import annotations

import argparse
import bisect
import sys

VALID_BASES = set("ACGT")


def load_gene_intervals(gff3_path: str) -> dict[str, list[tuple]]:
    """chrom -> sorted list of (start, end), 1-based inclusive, from `gene` features."""
    by_chrom: dict[str, list[tuple]] = {}
    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5 or fields[2] != "gene":
                continue
            by_chrom.setdefault(fields[0], []).append((int(fields[3]), int(fields[4])))
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    return by_chrom


def overlaps_gene(chrom: str, pos: int, gene_index: dict[str, list[tuple]]) -> bool:
    intervals = gene_index.get(chrom)
    if not intervals:
        return False
    starts = [iv[0] for iv in intervals]
    i = bisect.bisect_right(starts, pos) - 1
    while i >= 0 and intervals[i][0] <= pos:
        if intervals[i][0] <= pos <= intervals[i][1]:
            return True
        i -= 1
    return False


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--syri-out", required=True)
    ap.add_argument("--genes-gff3", default=None,
                    help="restrict to SNPs overlapping a gene body in this HA GFF3")
    ap.add_argument("--sample-name", default="HA_HB")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    gene_index = load_gene_intervals(args.genes_gff3) if args.genes_gff3 else None

    records = []
    n_snp_rows = 0
    n_kept_alleles = 0
    with open(args.syri_out) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            if fields[10] != "SNP":
                continue
            n_snp_rows += 1
            ref_chr, ref_pos, ref_allele, alt_allele = fields[0], fields[1], fields[3], fields[4]
            if ref_allele.upper() not in VALID_BASES or alt_allele.upper() not in VALID_BASES:
                continue
            pos = int(ref_pos)
            if gene_index is not None and not overlaps_gene(ref_chr, pos, gene_index):
                continue
            n_kept_alleles += 1
            records.append((ref_chr, pos, ref_allele.upper(), alt_allele.upper()))

    records.sort(key=lambda r: (r[0], r[1]))

    with open(args.out, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##source=symbiomics/syri_to_phased_vcf.py\n")
        out.write(
            "##INFO=<ID=SRC,Number=1,Type=String,Description="
            "\"HA-vs-HB assembly alignment (SyRI); phase known by construction, not statistically inferred\">\n"
        )
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + args.sample_name + "\n")
        for chrom, pos, ref, alt in records:
            out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tSRC=syri\tGT\t0|1\n")

    sys.stderr.write(
        f"INFO  ~ [syri] {n_snp_rows:,} SNP rows in {args.syri_out}, "
        f"{n_kept_alleles:,} kept" +
        (f" (gene-overlapping only)" if gene_index is not None else "") +
        f" -> {len(records):,} VCF records\n"
    )
    if not records:
        sys.stderr.write(
            "WARN  ~ [syri] 0 phased sites written -- phASER will have nothing to phase. "
            "Check syri.out actually contains SNP rows and, if --genes-gff3 was given, "
            "that its chromosome names match syri.out's.\n"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
