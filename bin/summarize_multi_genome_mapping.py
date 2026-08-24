#!/usr/bin/env python3
"""Compare per-sample HISAT2 mapping rates across multiple candidate genomes.

Built for the host/symbiont case this pipeline is named after: the same
transcriptome reads aligned once per candidate genome (see
subworkflows/local/multi_genome_mapping), each alignment producing a HISAT2
`--new-summary` file. This script turns that pile of per-(sample, genome)
summaries into one matrix and one best-genome call per sample -- the read
origin call a dual-organism or multi-candidate-symbiont dataset actually
needs, rather than N separate mapping-rate numbers nobody compares.

"Best genome" is a coverage argmax, not a taxonomic identification: a low
overall rate on every genome means the reads don't belong to any of the
candidates offered, and a thin margin between the top two means the call is
genuinely ambiguous (contamination, a close paralog, or two genomes that are
both partially right, e.g. an unremoved organelle sequence shared by host and
symbiont). Both are surfaced rather than papered over with a silent argmax.
"""
from __future__ import annotations

import argparse
import json
import re
import sys

OVERALL_RATE_RE = re.compile(r"Overall alignment rate:\s*([\d.]+)%")


def parse_hisat2_summary(path: str) -> float | None:
    """Return the overall alignment rate (0-100) from a hisat2 --new-summary file."""
    with open(path) as fh:
        text = fh.read()
    m = OVERALL_RATE_RE.search(text)
    return float(m.group(1)) if m else None


def read_manifest(path: str):
    rows = []
    with open(path) as fh:
        header = None
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                continue
            rows.append(dict(zip(header, fields)))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True,
                    help="TSV: sample_id, genome_id, summary_file")
    ap.add_argument("--out-matrix", required=True)
    ap.add_argument("--out-best", required=True)
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--ambiguous-margin", type=float, default=5.0,
                    help="flag a sample as ambiguous when the top two genomes' "
                         "alignment rates are within this many percentage points")
    args = ap.parse_args()

    rows = read_manifest(args.manifest)
    if not rows:
        sys.stderr.write("ERROR ~ [multi_genome_map] manifest has no rows\n")
        return 1

    # sample_id -> genome_id -> rate
    matrix: dict[str, dict[str, float]] = {}
    genomes: list[str] = []
    unparsed: list[str] = []

    for row in rows:
        sample, genome, summary = row["sample_id"], row["genome_id"], row["summary_file"]
        if genome not in genomes:
            genomes.append(genome)
        rate = parse_hisat2_summary(summary)
        if rate is None:
            unparsed.append(f"{sample}/{genome}")
            rate = float("nan")
        matrix.setdefault(sample, {})[genome] = rate

    samples = sorted(matrix)
    genomes = sorted(genomes)

    with open(args.out_matrix, "w") as fh:
        fh.write("sample_id\t" + "\t".join(genomes) + "\n")
        for sample in samples:
            row_vals = matrix[sample]
            fh.write(sample + "\t" + "\t".join(
                ("%.2f" % row_vals[g]) if g in row_vals and row_vals[g] == row_vals[g] else "NA"
                for g in genomes
            ) + "\n")

    best_rows = []
    with open(args.out_best, "w") as fh:
        fh.write("sample_id\tbest_genome\tbest_rate\trunner_up_genome\trunner_up_rate\tmargin\tambiguous\n")
        for sample in samples:
            ranked = sorted(
                ((g, r) for g, r in matrix[sample].items() if r == r),  # drop NaN
                key=lambda kv: kv[1], reverse=True,
            )
            if not ranked:
                fh.write(f"{sample}\tNA\tNA\tNA\tNA\tNA\tyes\n")
                best_rows.append({"sample_id": sample, "best_genome": None})
                continue
            best_genome, best_rate = ranked[0]
            if len(ranked) > 1:
                runner_genome, runner_rate = ranked[1]
                margin = best_rate - runner_rate
            else:
                runner_genome, runner_rate, margin = None, None, best_rate
            ambiguous = margin < args.ambiguous_margin
            fh.write(
                f"{sample}\t{best_genome}\t{best_rate:.2f}\t"
                f"{runner_genome or 'NA'}\t{('%.2f' % runner_rate) if runner_rate is not None else 'NA'}\t"
                f"{margin:.2f}\t{'yes' if ambiguous else 'no'}\n"
            )
            best_rows.append({
                "sample_id": sample,
                "best_genome": best_genome,
                "best_rate": best_rate,
                "runner_up_genome": runner_genome,
                "runner_up_rate": runner_rate,
                "margin": margin,
                "ambiguous": ambiguous,
            })
            if ambiguous:
                sys.stderr.write(
                    f"WARN  ~ [multi_genome_map] {sample}: {best_genome} "
                    f"({best_rate:.2f}%) leads {runner_genome} ({runner_rate:.2f}%) "
                    f"by only {margin:.2f} points -- ambiguous read origin\n"
                )

    summary = {
        "genomes": genomes,
        "samples": best_rows,
        "ambiguous_margin": args.ambiguous_margin,
        "unparsed_summaries": unparsed,
    }
    with open(args.out_json, "w") as fh:
        json.dump(summary, fh, indent=2)

    if unparsed:
        sys.stderr.write(
            f"WARN  ~ [multi_genome_map] could not parse an overall alignment "
            f"rate from {len(unparsed)} summary file(s): {', '.join(unparsed)}\n"
        )

    sys.stderr.write(
        f"INFO  ~ [multi_genome_map] {len(samples)} samples x {len(genomes)} "
        f"genomes -> {args.out_matrix}\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
