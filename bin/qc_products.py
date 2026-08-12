#!/usr/bin/env python3
"""Report the %hypothetical-protein metric from assign_products.py's output.

This is the QC guard that survives from the original AHRD-based design, minus
AHRD: a high fraction of proteins with no assigned product is a signal the
functional-annotation pass under-performed (wrong taxon, missing DB, too
divergent a Swiss-Prot reference), not something to bury in a log.
"""
from __future__ import annotations

import argparse
import json
import sys


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--products", required=True, help="assign_products.py output TSV")
    ap.add_argument("--max-hypothetical-frac", type=float, default=0.55)
    ap.add_argument("--out", default="products_qc.json")
    args = ap.parse_args()

    total = 0
    by_source: dict[str, int] = {}
    with open(args.products) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            total += 1
            source = fields[idx["source"]]
            by_source[source] = by_source.get(source, 0) + 1

    if total == 0:
        sys.stderr.write(f"ERROR ~ [products] {args.products} has no data rows\n")
        return 1

    hypothetical_frac = by_source.get("none", 0) / total
    result = {
        "total_proteins": total,
        "by_source": by_source,
        "hypothetical_fraction": round(hypothetical_frac, 4),
        "max_hypothetical_frac": args.max_hypothetical_frac,
        "over_threshold": hypothetical_frac > args.max_hypothetical_frac,
    }
    with open(args.out, "w") as fh:
        json.dump(result, fh, indent=2)

    sys.stderr.write(
        f"INFO  ~ [products] {hypothetical_frac*100:.1f}% hypothetical protein "
        f"({by_source.get('none', 0):,}/{total:,})\n"
    )
    if result["over_threshold"]:
        sys.stderr.write(
            f"WARN  ~ [products] hypothetical fraction {hypothetical_frac*100:.1f}% "
            f"exceeds the {args.max_hypothetical_frac*100:.0f}% guard. Check taxon, "
            f"reference proteome relevance, and whether Phase B modules "
            f"(InterProScan/KofamScan) would help -- see docs/functional_annotation.md.\n"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
