#!/usr/bin/env python3
"""Unit tests for bin/summarize_multi_genome_mapping.py. No pytest dependency:

    python3 tests/test_summarize_multi_genome_mapping.py
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(os.path.dirname(HERE), "bin", "summarize_multi_genome_mapping.py")

PASS, FAIL = 0, 0


def check(name: str, condition: bool, detail: str = "") -> None:
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  ok   {name}")
    else:
        FAIL += 1
        print(f"  FAIL {name}  {detail}")


def write_summary(tmp: str, name: str, rate: float) -> str:
    path = os.path.join(tmp, name)
    with open(path, "w") as fh:
        fh.write(
            "HISAT2 summary stats:\n"
            "\tTotal reads: 1000\n"
            "\t\tAligned 0 time: 10 (1.00%)\n"
            "\t\tAligned exactly 1 time: 900 (90.00%)\n"
            f"Overall alignment rate: {rate:.2f}%\n"
        )
    return path


def run(manifest: str, tmp: str, *extra: str):
    return subprocess.run(
        [sys.executable, SCRIPT, "--manifest", manifest,
         "--out-matrix", os.path.join(tmp, "matrix.tsv"),
         "--out-best", os.path.join(tmp, "best.tsv"),
         "--out-json", os.path.join(tmp, "summary.json"), *extra],
        capture_output=True, text=True,
    )


def write_manifest(tmp: str, rows: list[tuple[str, str, str]]) -> str:
    path = os.path.join(tmp, "manifest.tsv")
    with open(path, "w") as fh:
        fh.write("sample_id\tgenome_id\tsummary_file\n")
        for sample, genome, summary in rows:
            fh.write(f"{sample}\t{genome}\t{summary}\n")
    return path


def test_clear_winner_is_not_ambiguous():
    print("a sample that clearly favors one genome is not flagged ambiguous")
    with tempfile.TemporaryDirectory() as tmp:
        host = write_summary(tmp, "s1__host.summary.txt", 95.0)
        symb = write_summary(tmp, "s1__symbiont.summary.txt", 12.0)
        manifest = write_manifest(tmp, [
            ("s1", "host", host), ("s1", "symbiont", symb),
        ])
        proc = run(manifest, tmp)
        check("exit 0", proc.returncode == 0, proc.stderr)
        with open(os.path.join(tmp, "best.tsv")) as fh:
            rows = fh.readlines()
        fields = rows[1].rstrip("\n").split("\t")
        check("best genome is host", fields[1] == "host", fields)
        check("not ambiguous", fields[-1] == "no", fields)


def test_close_call_is_ambiguous():
    print("two genomes within the margin are flagged ambiguous")
    with tempfile.TemporaryDirectory() as tmp:
        a = write_summary(tmp, "s1__a.summary.txt", 60.0)
        b = write_summary(tmp, "s1__b.summary.txt", 58.0)
        manifest = write_manifest(tmp, [("s1", "a", a), ("s1", "b", b)])
        proc = run(manifest, tmp, "--ambiguous-margin", "5.0")
        check("exit 0", proc.returncode == 0, proc.stderr)
        with open(os.path.join(tmp, "best.tsv")) as fh:
            rows = fh.readlines()
        fields = rows[1].rstrip("\n").split("\t")
        check("ambiguous", fields[-1] == "yes", fields)


def test_matrix_has_one_row_per_sample_one_column_per_genome():
    print("matrix shape matches sample/genome counts, missing cells are NA")
    with tempfile.TemporaryDirectory() as tmp:
        s1a = write_summary(tmp, "s1__a.summary.txt", 80.0)
        s1b = write_summary(tmp, "s1__b.summary.txt", 20.0)
        s2a = write_summary(tmp, "s2__a.summary.txt", 30.0)
        manifest = write_manifest(tmp, [
            ("s1", "a", s1a), ("s1", "b", s1b), ("s2", "a", s2a),
        ])
        proc = run(manifest, tmp)
        check("exit 0", proc.returncode == 0, proc.stderr)
        with open(os.path.join(tmp, "matrix.tsv")) as fh:
            lines = [l.rstrip("\n").split("\t") for l in fh]
        check("header has both genomes", lines[0] == ["sample_id", "a", "b"], lines[0])
        s2_row = next(r for r in lines[1:] if r[0] == "s2")
        check("s2 missing genome b is NA", s2_row[2] == "NA", s2_row)


def test_unparsable_summary_reported_not_crashed():
    print("a summary file without an overall-rate line is reported, not fatal")
    with tempfile.TemporaryDirectory() as tmp:
        good = write_summary(tmp, "s1__a.summary.txt", 70.0)
        bad_path = os.path.join(tmp, "s1__b.summary.txt")
        with open(bad_path, "w") as fh:
            fh.write("not a real hisat2 summary\n")
        manifest = write_manifest(tmp, [("s1", "a", good), ("s1", "b", bad_path)])
        proc = run(manifest, tmp)
        check("exit 0", proc.returncode == 0, proc.stderr)
        check("unparsed noted in stderr", "could not parse" in proc.stderr, proc.stderr)
        with open(os.path.join(tmp, "summary.json")) as fh:
            payload = json.load(fh)
        check("unparsed_summaries recorded", len(payload["unparsed_summaries"]) == 1, payload)


def main() -> int:
    for test in [
        test_clear_winner_is_not_ambiguous,
        test_close_call_is_ambiguous,
        test_matrix_has_one_row_per_sample_one_column_per_genome,
        test_unparsable_summary_reported_not_crashed,
    ]:
        test()
    print(f"\n{PASS} passed, {FAIL} failed")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
