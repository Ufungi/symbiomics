#!/usr/bin/env python3
"""Split a combined HA+HB Salmon gene matrix (docs/haplotypes.md option 2)
into a summed total-expression matrix and a per-haplotype allele-specific
matrix, using the HA<->HB pairing table from syri_to_allele_table.py.

Input: salmon_to_gene_matrix.py's output run against the concatenated
HA.CDS.fa + HB.CDS.fa transcriptome (so gene_id values are a mix of HA and HB
gene IDs, whichever the CDS FASTA and GFF3 used).

Output:
  - summed.tsv: one row per allele PAIR, gene_id = "<ha_id>|<hb_id>", value =
    ha_count + hb_count. This is the PAV-inclusive, reference-bias-free total
    expression matrix -- the answer to "how much of this gene, from either
    haplotype".
  - allele_specific.tsv: unchanged, one row per original HA or HB gene_id, but
    now carrying a `haplotype` and `pair_id` column so a consumer can group by
    pair without re-joining the pairing table.
  - unpaired.tsv: genes that had no partner in the pairing table (an HA-only
    or HB-only gene -- a real PAV, or a synteny-block mismatch from
    syri_to_allele_table.py; the two are not distinguished here without the
    table's own block-mismatch report).
"""
from __future__ import annotations

import argparse
import sys


def read_matrix(path: str) -> tuple[list[str], dict[str, list[str]]]:
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        samples = header[1:]
        rows = {}
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            rows[fields[0]] = fields[1:]
    return samples, rows


def read_pairs(path: str) -> tuple[dict[str, str], dict[str, str]]:
    """Returns (ha_to_hb, hb_to_ha), skipping rows where either side is blank."""
    ha_to_hb, hb_to_ha = {}, {}
    with open(path) as fh:
        header = fh.readline()
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            ha_id, hb_id = fields[0], fields[1]
            if ha_id and hb_id:
                ha_to_hb[ha_id] = hb_id
                hb_to_ha[hb_id] = ha_id
    return ha_to_hb, hb_to_ha


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix", required=True, help="combined salmon gene matrix (counts or TPM)")
    ap.add_argument("--pairs", required=True, help="syri_to_allele_table.py output")
    ap.add_argument("--out-summed", required=True)
    ap.add_argument("--out-allele-specific", required=True)
    ap.add_argument("--out-unpaired", required=True)
    args = ap.parse_args()

    samples, rows = read_matrix(args.matrix)
    ha_to_hb, hb_to_ha = read_pairs(args.pairs)

    # A gene named in the pairing table but absent from the matrix is not
    # necessarily unquantified: Salmon deduplicates byte-identical reference
    # sequences at index time and reports only one representative, so a very
    # highly conserved allele pair can legitimately have only one side survive
    # in quant.sf. Treating the missing side as zero is correct either way --
    # if it was deduplicated, the surviving id's count already carries the
    # combined signal (Salmon could not have told the reads apart either), so
    # summing it against zero reproduces that combined figure exactly; if it
    # was genuinely unquantified, zero is simply the right answer.
    all_ids = set(rows)
    paired_ha = set(ha_to_hb) & all_ids
    paired_hb = set(hb_to_ha) & all_ids
    n_samples = len(samples)
    zero_row = ["0"] * n_samples

    seen_pairs = set()
    n_deduped = 0
    with open(args.out_summed, "w") as summed, \
         open(args.out_allele_specific, "w") as allele, \
         open(args.out_unpaired, "w") as unpaired:

        summed.write("pair_id\t" + "\t".join(samples) + "\n")
        allele.write("gene_id\thaplotype\tpair_id\t" + "\t".join(samples) + "\n")
        unpaired.write("gene_id\thaplotype\t" + "\t".join(samples) + "\n")

        for ha_id in sorted(paired_ha):
            hb_id = ha_to_hb[ha_id]
            pair_id = f"{ha_id}|{hb_id}"
            if pair_id in seen_pairs:
                continue
            seen_pairs.add(pair_id)
            if hb_id not in rows:
                n_deduped += 1
            ha_row = rows[ha_id]
            hb_row = rows.get(hb_id, zero_row)
            ha_vals = [float(v) for v in ha_row]
            hb_vals = [float(v) for v in hb_row]
            summed.write(pair_id + "\t" + "\t".join(f"{a+b:.3f}" for a, b in zip(ha_vals, hb_vals)) + "\n")
            allele.write(f"{ha_id}\tHA\t{pair_id}\t" + "\t".join(ha_row) + "\n")
            allele.write(f"{hb_id}\tHB\t{pair_id}\t" + "\t".join(hb_row) + "\n")

        for gene_id in sorted(all_ids - paired_ha - paired_hb):
            hap = "HA" if gene_id in ha_to_hb or gene_id not in hb_to_ha else "HB"
            unpaired.write(f"{gene_id}\t{hap}\t" + "\t".join(rows[gene_id]) + "\n")

    n_pairs = len(seen_pairs)
    n_unpaired = len(all_ids - paired_ha - paired_hb)
    sys.stderr.write(
        f"INFO  ~ [haplotype-split] {n_pairs:,} allele pairs summed"
        + (f" ({n_deduped:,} with a partner Salmon deduplicated away -- treated as zero, "
           f"see the script docstring)" if n_deduped else "")
        + f", {n_unpaired:,} genes unpaired (of {len(all_ids):,} total gene_ids in the matrix)\n"
    )
    if len(all_ids) and n_unpaired / len(all_ids) > 0.3:
        sys.stderr.write(
            f"WARN  ~ [haplotype-split] {n_unpaired/len(all_ids)*100:.0f}% of genes are "
            f"unpaired -- check that --pairs matches this matrix's gene_id convention "
            f"(same GFF3s used for both), or that the SyRI alignment covers most of the genome.\n"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
