#!/usr/bin/env python3
"""Build ONE diploid CDS FASTA (HA+HB) with SEESAW-compatible allele suffixes.

SEESAW allelic-imbalance analysis needs the two haplotypes of a diploid
assembly quantified together but kept apart: every transcript must carry a
consistent allele suffix so the two copies of a gene can be told apart, and it
needs a `transcript_id -> gene_id` table to aggregate transcript-level
quantification back to genes. The two P. densiflora haplotype CDS files name
their genes identically except for a trailing haplotype letter (Pd01G00010A vs
Pd01G00010B), which no downstream tool understands as a pair, so here we:

  1. rename every HA transcript  PdNNNgNNNNNA -> PdNNNgNNNNN_A
     and every HB transcript  PdNNNgNNNNNB -> PdNNNgNNNNN_B,
  2. concatenate the two renamed sets into a single CDS FASTA,
  3. emit a tx2gene table (identity -- one CDS per gene, keyed by gene ID; the
     Pinus convention the pipeline's GFF3s already follow), and
  4. emit a pairing table (stem -> ha_id, hb_id) for the genes present in both
     haplotypes, so a consumer can sum or compare the two alleles.

The rename step is exactly the hazard `check_disjoint_fasta_ids.py` guards
against for arbitrary concatenations: two files sharing an ID would silently
merge distinct sequences under one name. Here the collision is deliberate and
checked -- the suffixes must make the combined ID space *disjoint*, so we verify
that no renamed ID repeats within or across the two files and refuse rather
than guess if it does (see docs/haplotypes.md).
"""
from __future__ import annotations

import argparse
import sys


def scan(path: str, allele: str) -> tuple[list[str], set[str], list[str]]:
    """Return (renamed_ids, stems, offenders) for one haplotype CDS FASTA.

    `allele` is the single trailing letter already present on every ID ("A" or
    "B"). The rename turns that suffix into an underscore form: id[:-1] + "_" +
    allele. IDs not ending in `allele` are collected as offenders (not renamed).
    """
    renamed: list[str] = []
    stems: set[str] = set()
    offenders: list[str] = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            seq_id = line[1:].split()[0]
            if not seq_id.endswith(allele):
                offenders.append(seq_id)
                continue
            stem = seq_id[:-1]
            stems.add(stem)
            renamed.append(f"{stem}_{allele}")
    return renamed, stems, offenders


def write_renamed(out_fh, path: str, allele: str) -> None:
    """Stream one FASTA, rewriting header IDs and passing sequences through."""
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                fields = line[1:].rstrip("\n").split(None, 1)
                rest = f" {fields[1]}" if len(fields) > 1 else ""
                out_fh.write(f">{fields[0][:-1]}_{allele}{rest}\n")
            else:
                out_fh.write(line)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ha", required=True, help="haplotype A CDS FASTA")
    ap.add_argument("--hb", required=True, help="haplotype B CDS FASTA")
    ap.add_argument("--out-fasta", required=True, help="diploid CDS FASTA (_A + _B)")
    ap.add_argument("--out-tx2gene", required=True, help="transcript_id -> gene_id TSV")
    ap.add_argument("--out-pairs", required=True, help="stem -> ha_id, hb_id TSV")
    args = ap.parse_args()

    ha_ids, ha_stems, ha_bad = scan(args.ha, "A")
    hb_ids, hb_stems, hb_bad = scan(args.hb, "B")

    if ha_bad or hb_bad:
        for label, path, bad in (("HA", args.ha, ha_bad), ("HB", args.hb, hb_bad)):
            if bad:
                sys.stderr.write(
                    f"ERROR ~ [prepare-diploid] {len(bad)} {label} ID(s) in {path} "
                    f"do not end with the expected haplotype letter '{label}':\n"
                )
                for seq_id in bad[:10]:
                    sys.stderr.write(f"        {seq_id!r}\n")
                if len(bad) > 10:
                    sys.stderr.write(f"        ... and {len(bad) - 10} more\n")
        sys.stderr.write(
            "        The rename assumes ALL HA IDs end in 'A' and ALL HB IDs in "
            "'B'; a mismatch would produce an inconsistent _A/_B namespace.\n"
        )
        return 1

    # The rename must make the combined ID space disjoint -- the same hazard
    # check_disjoint_fasta_ids.py refuses on. Collisions within a file (a
    # repeated stem) or across files are both fatal.
    seen: set[str] = set()
    collisions: list[str] = []
    for new_id in ha_ids + hb_ids:
        if new_id in seen:
            collisions.append(new_id)
        else:
            seen.add(new_id)
    if collisions:
        sys.stderr.write(
            f"ERROR ~ [prepare-diploid] {len(collisions)} renamed ID(s) collide "
            f"within or across the two haplotype FASTAs:\n"
        )
        for new_id in collisions[:10]:
            sys.stderr.write(f"        {new_id!r}\n")
        if len(collisions) > 10:
            sys.stderr.write(f"        ... and {len(collisions) - 10} more\n")
        sys.stderr.write(
            "        Concatenating these would silently merge two different "
            "sequences under one name (see docs/haplotypes.md).\n"
        )
        return 1

    with open(args.out_fasta, "w", encoding="utf-8", newline="\n") as out_fa, \
         open(args.out_tx2gene, "w", encoding="utf-8", newline="\n") as out_tx, \
         open(args.out_pairs, "w", encoding="utf-8", newline="\n") as out_pairs:
        write_renamed(out_fa, args.ha, "A")
        write_renamed(out_fa, args.hb, "B")

        out_tx.write("transcript_id\tgene_id\n")
        for new_id in ha_ids + hb_ids:
            out_tx.write(f"{new_id}\t{new_id}\n")

        paired = ha_stems & hb_stems
        out_pairs.write("stem\tha_id\thb_id\n")
        for stem in sorted(paired):
            out_pairs.write(f"{stem}\t{stem}_A\t{stem}_B\n")

    ha_only = len(ha_stems - hb_stems)
    hb_only = len(hb_stems - ha_stems)
    sys.stderr.write(
        f"INFO  ~ [prepare-diploid] HA {len(ha_ids):,} + HB {len(hb_ids):,} = "
        f"{len(ha_ids) + len(hb_ids):,} transcripts; {len(paired):,} paired stems, "
        f"{ha_only:,} HA-only, {hb_only:,} HB-only -> {args.out_fasta}\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
