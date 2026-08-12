#!/usr/bin/env python3
"""Build an HA<->HB gene-level allele-pairing table from SyRI's syntenic
(SYN) blocks.

For each SYN block, SyRI gives a co-linear correspondence between a ref (HA)
interval and a query (HB) interval. Within one block, genes are expected to
appear in the same relative order on both sides (that is what "syntenic"
means) -- so genes are paired positionally within each block: the Nth HA gene
overlapping the block pairs with the Nth HB gene overlapping the block.

This is a standard, defensible approximation for allele pairing between two
haplotypes of the same individual (not two different species, where local
rearrangements would break the assumption more often). It is NOT attempted
across blocks, and a block whose HA/HB gene counts disagree reports its genes
as unpaired rather than guessing which one is missing.

Used by:
  - the diploid Salmon EM split (docs/haplotypes.md option 2)
  - nothing in the phASER route directly (that uses the SNP-level VCF from
    syri_to_phased_vcf.py) -- this table is for gene-level bookkeeping only.
"""
from __future__ import annotations

import argparse
import sys


def load_genes(gff3_path: str) -> list[tuple]:
    """Returns [(chrom, start, end, gene_id), ...] from `gene` features."""
    genes = []
    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            gid = None
            for kv in fields[8].split(";"):
                kv = kv.strip()
                if kv.startswith("ID="):
                    gid = kv[3:]
                    break
            if gid:
                genes.append((fields[0], int(fields[3]), int(fields[4]), gid))
    genes.sort()
    return genes


def load_syn_blocks(syri_out: str) -> list[tuple]:
    """Returns [(ref_chr, ref_start, ref_end, query_chr, query_start, query_end), ...]."""
    blocks = []
    with open(syri_out) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11 or fields[10] != "SYN":
                continue
            try:
                blocks.append((
                    fields[0], int(fields[1]), int(fields[2]),
                    fields[5], int(fields[6]), int(fields[7]),
                ))
            except ValueError:
                continue
    return blocks


def genes_in_interval(genes: list[tuple], chrom: str, start: int, end: int) -> list[tuple]:
    lo, hi = min(start, end), max(start, end)
    return [g for g in genes if g[0] == chrom and g[1] >= lo and g[2] <= hi]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--syri-out", required=True)
    ap.add_argument("--ha-gff3", required=True)
    ap.add_argument("--hb-gff3", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    ha_genes = load_genes(args.ha_gff3)
    hb_genes = load_genes(args.hb_gff3)
    blocks = load_syn_blocks(args.syri_out)
    sys.stderr.write(
        f"INFO  ~ [pairing] {len(ha_genes):,} HA genes, {len(hb_genes):,} HB genes, "
        f"{len(blocks):,} syntenic blocks\n"
    )

    paired = 0
    unpaired_mismatch = 0
    ha_seen: set[str] = set()
    hb_seen: set[str] = set()

    with open(args.out, "w") as out:
        out.write("ha_gene_id\thb_gene_id\tblock_ref\tblock_query\n")
        for ref_chr, ref_start, ref_end, query_chr, query_start, query_end in blocks:
            ha_in_block = genes_in_interval(ha_genes, ref_chr, ref_start, ref_end)
            hb_in_block = genes_in_interval(hb_genes, query_chr, query_start, query_end)
            reversed_block = query_start > query_end
            if reversed_block:
                hb_in_block = sorted(hb_in_block, key=lambda g: -g[1])
            block_label = f"{ref_chr}:{ref_start}-{ref_end}"
            block_label_q = f"{query_chr}:{min(query_start,query_end)}-{max(query_start,query_end)}"

            if len(ha_in_block) == len(hb_in_block) and ha_in_block:
                for (_, _, _, ha_id), (_, _, _, hb_id) in zip(ha_in_block, hb_in_block):
                    out.write(f"{ha_id}\t{hb_id}\t{block_label}\t{block_label_q}\n")
                    ha_seen.add(ha_id)
                    hb_seen.add(hb_id)
                    paired += 1
            else:
                unpaired_mismatch += len(ha_in_block) + len(hb_in_block)

        ha_unpaired = [g[3] for g in ha_genes if g[3] not in ha_seen]
        for gid in ha_unpaired:
            out.write(f"{gid}\t\t\t\n")
        hb_unpaired = [g[3] for g in hb_genes if g[3] not in hb_seen]
        for gid in hb_unpaired:
            out.write(f"\t{gid}\t\t\n")

    sys.stderr.write(
        f"INFO  ~ [pairing] {paired:,} gene pairs from collinear blocks; "
        f"{len(ha_unpaired):,} HA and {len(hb_unpaired):,} HB genes unpaired "
        f"({unpaired_mismatch:,} genes sat in blocks with mismatched HA/HB counts) "
        f"-> {args.out}\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
