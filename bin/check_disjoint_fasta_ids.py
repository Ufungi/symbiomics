#!/usr/bin/env python3
"""Concatenate FASTA files after confirming their sequence IDs are disjoint.

This is the CDS-level version of a trap already documented for genomes in
docs/haplotypes.md: if two haplotypes' CDS files reuse an ID (or, worse, name a
gene identically because their source GFF3s used the same numbering scheme),
concatenating them silently merges two different sequences under one name and
every downstream count becomes meaningless. Refuse rather than guess.
"""
from __future__ import annotations

import argparse
import sys


def read_ids(path: str) -> set[str]:
    ids = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.add(line[1:].split()[0])
    return ids


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("fastas", nargs="+")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if len(args.fastas) < 2:
        # Nothing to check; just copy through.
        with open(args.out, "w") as out:
            with open(args.fastas[0]) as fh:
                out.write(fh.read())
        return 0

    seen: dict[str, str] = {}
    collisions = []
    for path in args.fastas:
        ids = read_ids(path)
        for seq_id in ids:
            if seq_id in seen:
                collisions.append((seq_id, seen[seq_id], path))
            else:
                seen[seq_id] = path

    if collisions:
        sys.stderr.write(
            f"ERROR ~ [salmon] {len(collisions)} sequence ID(s) are shared across "
            f"the inputs meant to represent different haplotypes/sources:\n"
        )
        for seq_id, first, second in collisions[:10]:
            sys.stderr.write(f"        {seq_id!r}: in both {first} and {second}\n")
        if len(collisions) > 10:
            sys.stderr.write(f"        ... and {len(collisions) - 10} more\n")
        sys.stderr.write(
            "        Concatenating these would silently merge two different "
            "sequences under one name (the same hazard as two haplotypes both "
            "naming a chromosome 'chr1a' -- see docs/haplotypes.md).\n"
            "        Rename the IDs to be unique per source before combining them.\n"
        )
        return 1

    with open(args.out, "w") as out:
        for path in args.fastas:
            with open(path) as fh:
                out.write(fh.read())

    sys.stderr.write(
        f"INFO  ~ [salmon] {len(seen):,} disjoint sequence IDs across "
        f"{len(args.fastas)} file(s) -> {args.out}\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
