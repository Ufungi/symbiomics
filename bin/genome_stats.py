#!/usr/bin/env python3
"""Single-pass genome statistics for the policy engine.

Streams the FASTA once and reports total bp, contig count, longest contig, N
fraction and softmask (lowercase) fraction. Softmask fraction is what tells the
pipeline whether --premasked is honest: a genome that claims to be masked but
has 0% lowercase is not masked.

Emits JSON so decide_strategy.py can consume it without re-parsing FASTA.
"""
from __future__ import annotations

import argparse
import json
import sys

LOWER = set(b"acgtn")
UPPER = set(b"ACGTN")


def scan(path: str, sample_bp: int = 0):
    total = 0
    contigs = 0
    longest = 0
    n_count = 0
    lower = 0
    current = 0
    scanned = 0
    truncated = False

    with open(path, "rb") as fh:
        for raw in fh:
            if raw.startswith(b">"):
                contigs += 1
                if current > longest:
                    longest = current
                total += current
                current = 0
                continue
            line = raw.rstrip()
            current += len(line)
            # Composition is measured on a prefix when sample_bp is set; length
            # accounting always covers the whole file.
            if not truncated:
                scanned += len(line)
                n_count += line.count(b"N") + line.count(b"n")
                lower += sum(1 for b in line if 97 <= b <= 122)
                if sample_bp and scanned >= sample_bp:
                    truncated = True

    if current > longest:
        longest = current
    total += current

    return {
        "total_bp": total,
        "contigs": contigs,
        "max_contig_bp": longest,
        "n_fraction": round(n_count / scanned, 6) if scanned else 0.0,
        "softmask_fraction": round(lower / scanned, 6) if scanned else 0.0,
        "composition_scanned_bp": scanned,
        "composition_is_sampled": truncated,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", default="genome_stats.json")
    ap.add_argument(
        "--sample-bp", type=int, default=200_000_000,
        help="measure N/softmask fraction on this many bp (0 = whole genome). "
             "Length accounting is always exact.",
    )
    args = ap.parse_args()

    stats = scan(args.fasta, args.sample_bp)
    if stats["total_bp"] == 0:
        sys.stderr.write(f"ERROR ~ [genome] {args.fasta} contains no sequence\n")
        return 1

    with open(args.out, "w") as fh:
        json.dump(stats, fh, indent=2)

    gb = stats["total_bp"] / 1e9
    sys.stderr.write(
        f"INFO  ~ [genome] {stats['total_bp']:,} bp ({gb:.3f} Gb), "
        f"{stats['contigs']:,} contigs, longest {stats['max_contig_bp']:,} bp, "
        f"softmasked {stats['softmask_fraction'] * 100:.1f}%\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
