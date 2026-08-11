#!/usr/bin/env python3
"""Work out which feature type and grouping attribute to count on.

Published GFF3 files vary far more than counting tools assume. featureCounts
defaults to `-t exon -g gene_id`, and a file that has no `gene_id` attribute
fails with an error that says nothing useful about the cause.

A real example -- the published Pinus densiflora v1.0 annotation -- carries only
`ID=` on every line, no `Parent=`, no `gene_id`, and repeats the gene's own ID on
its exon and CDS children. `-t exon -g ID` counts it correctly; the default does
not run at all.

This inspects the file and picks a combination that will work, then reports what
it chose and why. It never guesses silently.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import json
import re
import sys

# Preference order for the grouping attribute. gene_id first because it is what
# every downstream convention expects when it is present.
ATTR_PREFERENCE = ["gene_id", "gene", "ID", "Name", "locus_tag", "Parent"]
FEATURE_PREFERENCE = ["exon", "CDS", "cds", "transcript", "mRNA", "gene"]


def opener(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_attrs(field: str) -> dict:
    out = {}
    # GTF: key "value";   GFF3: key=value;
    if "=" in field and not re.search(r'\w\s+"', field):
        for kv in field.split(";"):
            kv = kv.strip()
            if "=" in kv:
                key, _, value = kv.partition("=")
                out[key.strip()] = value.strip()
    else:
        for kv in field.split(";"):
            kv = kv.strip()
            match = re.match(r'^(\S+)\s+"?([^"]*)"?$', kv)
            if match:
                out[match.group(1)] = match.group(2)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--out", default="annotation_spec.json")
    ap.add_argument("--want-feature", default="exon")
    ap.add_argument("--want-attribute", default="gene_id")
    ap.add_argument("--max-lines", type=int, default=500000)
    args = ap.parse_args()

    feature_counts = collections.Counter()
    attrs_by_feature = collections.defaultdict(collections.Counter)
    values_by_pair = collections.defaultdict(set)

    with opener(args.annotation) as fh:
        for i, line in enumerate(fh):
            if i >= args.max_lines:
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            feature_counts[feature] += 1
            attrs = parse_attrs(fields[8])
            for key, value in attrs.items():
                attrs_by_feature[feature][key] += 1
                if len(values_by_pair[(feature, key)]) < 200000:
                    values_by_pair[(feature, key)].add(value)

    if not feature_counts:
        sys.stderr.write(
            f"ERROR ~ [annotation] {args.annotation} has no usable feature lines.\n"
        )
        return 1

    notes = []

    # --- feature type -------------------------------------------------------
    feature = args.want_feature
    if feature_counts.get(feature, 0) == 0:
        for candidate in FEATURE_PREFERENCE:
            if feature_counts.get(candidate, 0) > 0:
                notes.append(
                    f"requested feature '{args.want_feature}' is absent; "
                    f"using '{candidate}' ({feature_counts[candidate]:,} lines)"
                )
                feature = candidate
                break
        else:
            feature = feature_counts.most_common(1)[0][0]
            notes.append(f"falling back to the commonest feature '{feature}'")

    # Case-variant features (a real file had both CDS and cds) would be counted
    # only in part; say so rather than silently dropping the minority spelling.
    variants = [f for f in feature_counts
                if f.lower() == feature.lower() and f != feature]
    if variants:
        notes.append(
            f"the file also uses {variants} -- case variants of '{feature}'; "
            f"{sum(feature_counts[v] for v in variants):,} lines are not counted"
        )

    # --- grouping attribute -------------------------------------------------
    present = attrs_by_feature[feature]
    n_feature = feature_counts[feature]
    attribute = None
    if present.get(args.want_attribute, 0) >= n_feature * 0.99:
        attribute = args.want_attribute
    else:
        if args.want_attribute in present:
            notes.append(
                f"attribute '{args.want_attribute}' is on only "
                f"{present[args.want_attribute]:,}/{n_feature:,} '{feature}' lines"
            )
        else:
            notes.append(
                f"attribute '{args.want_attribute}' is absent from '{feature}' "
                f"lines (present: {sorted(present)})"
            )
        for candidate in ATTR_PREFERENCE:
            if present.get(candidate, 0) >= n_feature * 0.99:
                attribute = candidate
                notes.append(f"using '{candidate}' instead")
                break

    if attribute is None:
        sys.stderr.write(
            f"ERROR ~ [annotation] no attribute covers every '{feature}' line in "
            f"{args.annotation}.\n"
            f"        attributes seen on '{feature}': "
            f"{dict(present.most_common(6))}\n"
            f"        Set --count_feature / --count_attribute explicitly.\n"
        )
        return 1

    n_groups = len(values_by_pair[(feature, attribute)])
    n_genes = feature_counts.get("gene", 0)
    if n_genes and abs(n_groups - n_genes) / n_genes > 0.05:
        notes.append(
            f"grouping by '{attribute}' yields {n_groups:,} groups but the file "
            f"declares {n_genes:,} genes -- check this is the level you want"
        )

    spec = {
        "annotation": args.annotation,
        "feature": feature,
        "attribute": attribute,
        "n_feature_lines": n_feature,
        "n_groups": n_groups,
        "n_genes_declared": n_genes,
        "feature_types": dict(feature_counts.most_common(10)),
        "notes": notes,
    }
    with open(args.out, "w") as fh:
        json.dump(spec, fh, indent=2)

    sys.stderr.write(
        f"INFO  ~ [annotation] counting '{feature}' grouped by '{attribute}' "
        f"-> {n_groups:,} features\n"
    )
    for note in notes:
        sys.stderr.write(f"WARN  ~ [annotation] {note}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
