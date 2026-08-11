#!/usr/bin/env python3
"""Unit tests for bin/infer_layout.py. No pytest dependency -- run directly:

    python3 tests/test_infer_layout.py
"""
from __future__ import annotations

import gzip
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(os.path.dirname(HERE), "bin", "infer_layout.py")

PASS, FAIL = 0, 0


def make_fastq(path: str, read_id: str = "READ", n: int = 4, mate: str = "1") -> None:
    with gzip.open(path, "wt") as fh:
        for i in range(n):
            fh.write(f"@{read_id}{i}/{mate}\nACGTACGTAC\n+\nIIIIIIIIII\n")


def run(sheet: str, tmp: str, *extra: str):
    return subprocess.run(
        [sys.executable, SCRIPT, "--input", sheet,
         "--out-tsv", os.path.join(tmp, "out.tsv"),
         "--out-json", os.path.join(tmp, "out.json"), *extra],
        capture_output=True, text=True,
    )


def write_sheet(tmp: str, body: str, name: str = "sheet.tsv") -> str:
    path = os.path.join(tmp, name)
    with open(path, "w") as fh:
        fh.write(body)
    return path


def check(name: str, condition: bool, detail: str = "") -> None:
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  ok   {name}")
    else:
        FAIL += 1
        print(f"  FAIL {name}  {detail}")


def resolved(tmp: str):
    rows = {}
    with open(os.path.join(tmp, "out.tsv")) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            values = line.rstrip("\n").split("\t")
            rows[values[0]] = dict(zip(header, values))
    return rows


def test_explicit_and_inference():
    print("explicit layout + filesystem inference")
    with tempfile.TemporaryDirectory() as tmp:
        make_fastq(os.path.join(tmp, "A_R1.fastq.gz"), "A", mate="1")
        make_fastq(os.path.join(tmp, "A_R2.fastq.gz"), "A", mate="2")
        make_fastq(os.path.join(tmp, "B.fastq.gz"), "B", mate="1")
        make_fastq(os.path.join(tmp, "C_1.fq.gz"), "C", mate="1")
        make_fastq(os.path.join(tmp, "C_2.fq.gz"), "C", mate="2")
        sheet = write_sheet(tmp, (
            "sample_id\tlayout\tfastq_1\tfastq_2\n"
            f"A\t-\t{tmp}/A_R1.fastq.gz\t-\n"
            f"B\t-\t{tmp}/B.fastq.gz\t-\n"
            f"C\tPE\t{tmp}/C_1.fq.gz\t{tmp}/C_2.fq.gz\n"
        ))
        proc = run(sheet, tmp)
        check("exit 0", proc.returncode == 0, proc.stderr)
        rows = resolved(tmp)
        check("A inferred PE via _R1->_R2", rows["A"]["layout"] == "PE", rows["A"])
        check("A mate path filled", rows["A"]["fastq_2"].endswith("A_R2.fastq.gz"))
        check("B stays SE", rows["B"]["layout"] == "SE", rows["B"])
        check("C explicit PE", rows["C"]["layout"] == "PE")


def test_se_with_mate_is_an_error():
    print("layout=SE while fastq_2 is set must fail")
    with tempfile.TemporaryDirectory() as tmp:
        make_fastq(os.path.join(tmp, "X_R1.fastq.gz"), "X", mate="1")
        make_fastq(os.path.join(tmp, "X_R2.fastq.gz"), "X", mate="2")
        sheet = write_sheet(tmp, (
            "sample_id\tlayout\tfastq_1\tfastq_2\n"
            f"X\tSE\t{tmp}/X_R1.fastq.gz\t{tmp}/X_R2.fastq.gz\n"
        ))
        proc = run(sheet, tmp)
        check("exit non-zero", proc.returncode != 0)
        check("message names the hazard",
              "Refusing to silently discard a mate file" in proc.stderr, proc.stderr)


def test_pe_without_mate_is_an_error():
    print("layout=PE without fastq_2 must fail")
    with tempfile.TemporaryDirectory() as tmp:
        make_fastq(os.path.join(tmp, "Y.fastq.gz"), "Y")
        sheet = write_sheet(tmp, (
            "sample_id\tlayout\tfastq_1\tfastq_2\n"
            f"Y\tPE\t{tmp}/Y.fastq.gz\t-\n"
        ))
        proc = run(sheet, tmp)
        check("exit non-zero", proc.returncode != 0)
        check("message explains", "layout=PE but fastq_2 is '-'" in proc.stderr, proc.stderr)


def test_mismatched_read_ids_rejected():
    print("inferred mate with disagreeing read IDs must fail")
    with tempfile.TemporaryDirectory() as tmp:
        make_fastq(os.path.join(tmp, "Z_R1.fastq.gz"), "AAA", mate="1")
        make_fastq(os.path.join(tmp, "Z_R2.fastq.gz"), "BBB", mate="2")
        sheet = write_sheet(tmp, (
            "sample_id\tlayout\tfastq_1\tfastq_2\n"
            f"Z\tauto\t{tmp}/Z_R1.fastq.gz\t-\n"
        ))
        proc = run(sheet, tmp)
        check("exit non-zero", proc.returncode != 0)
        check("message names read IDs", "read IDs disagree" in proc.stderr, proc.stderr)


def test_duplicate_sample_ids():
    print("duplicate sample_id must fail")
    with tempfile.TemporaryDirectory() as tmp:
        make_fastq(os.path.join(tmp, "D.fastq.gz"), "D")
        sheet = write_sheet(tmp, (
            "sample_id\tfastq_1\n"
            f"D\t{tmp}/D.fastq.gz\n"
            f"D\t{tmp}/D.fastq.gz\n"
        ))
        proc = run(sheet, tmp)
        check("exit non-zero", proc.returncode != 0)
        check("message cites first line", "is duplicated" in proc.stderr, proc.stderr)


def test_missing_file_message():
    print("missing fastq reports the resolved path")
    with tempfile.TemporaryDirectory() as tmp:
        sheet = write_sheet(tmp, "sample_id\tfastq_1\nE\tnope.fastq.gz\n")
        proc = run(sheet, tmp, "--raw-dir", tmp)
        check("exit non-zero", proc.returncode != 0)
        check("shows attempted path", f"{tmp}/nope.fastq.gz" in proc.stderr, proc.stderr)
        check("mentions raw_dir", "--raw_dir" in proc.stderr, proc.stderr)


def test_comments_aliases_and_replicates():
    print("comments, header aliases, technical replicates")
    with tempfile.TemporaryDirectory() as tmp:
        for tag in ("L1", "L2"):
            make_fastq(os.path.join(tmp, f"{tag}_R1.fastq.gz"), "M", mate="1")
            make_fastq(os.path.join(tmp, f"{tag}_R2.fastq.gz"), "M", mate="2")
        sheet = write_sheet(tmp, (
            "# a comment line\n"
            "\n"
            "sample\tR1\tR2\n"
            f"M\t{tmp}/L1_R1.fastq.gz,{tmp}/L2_R1.fastq.gz"
            f"\t{tmp}/L1_R2.fastq.gz,{tmp}/L2_R2.fastq.gz\n"
        ))
        proc = run(sheet, tmp)
        check("exit 0", proc.returncode == 0, proc.stderr)
        rows = resolved(tmp)
        check("aliases accepted", "M" in rows, rows)
        check("PE with 2 replicates", rows["M"]["fastq_1"].count(",") == 1, rows["M"])


def test_unequal_replicate_lists():
    print("unequal technical-replicate lists must fail")
    with tempfile.TemporaryDirectory() as tmp:
        for tag in ("L1", "L2"):
            make_fastq(os.path.join(tmp, f"{tag}_R1.fastq.gz"), "N", mate="1")
        make_fastq(os.path.join(tmp, "L1_R2.fastq.gz"), "N", mate="2")
        sheet = write_sheet(tmp, (
            "sample_id\tfastq_1\tfastq_2\n"
            f"N\t{tmp}/L1_R1.fastq.gz,{tmp}/L2_R1.fastq.gz\t{tmp}/L1_R2.fastq.gz\n"
        ))
        proc = run(sheet, tmp)
        check("exit non-zero", proc.returncode != 0)
        check("explains lengths", "same length" in proc.stderr, proc.stderr)


def test_bam_and_fastq_mutually_exclusive():
    print("bam + fastq_1 together must fail")
    with tempfile.TemporaryDirectory() as tmp:
        make_fastq(os.path.join(tmp, "P.fastq.gz"), "P")
        open(os.path.join(tmp, "P.bam"), "wb").close()
        sheet = write_sheet(tmp, (
            "sample_id\tfastq_1\tbam\n"
            f"P\t{tmp}/P.fastq.gz\t{tmp}/P.bam\n"
        ))
        proc = run(sheet, tmp)
        check("exit non-zero", proc.returncode != 0)
        check("explains exclusivity", "mutually exclusive" in proc.stderr, proc.stderr)


def test_csv_uses_semicolon_for_replicates():
    print("CSV input: ';' separates technical replicates")
    with tempfile.TemporaryDirectory() as tmp:
        make_fastq(os.path.join(tmp, "Q1_R1.fastq.gz"), "Q", mate="1")
        make_fastq(os.path.join(tmp, "Q2_R1.fastq.gz"), "Q", mate="1")
        sheet = write_sheet(tmp, (
            "sample_id,layout,fastq_1,fastq_2\n"
            f"Q,SE,{tmp}/Q1_R1.fastq.gz;{tmp}/Q2_R1.fastq.gz,-\n"
        ), name="sheet.csv")
        proc = run(sheet, tmp)
        check("exit 0", proc.returncode == 0, proc.stderr)
        rows = resolved(tmp)
        check("two replicate files", rows["Q"]["fastq_1"].count(",") == 1, rows["Q"])


def test_use_for_and_strandedness_validation():
    print("use_for / strandedness validation")
    with tempfile.TemporaryDirectory() as tmp:
        make_fastq(os.path.join(tmp, "R.fastq.gz"), "R")
        good = write_sheet(tmp, (
            "sample_id\tfastq_1\tstrandedness\tuse_for\n"
            f"R\t{tmp}/R.fastq.gz\treverse\tstringtie,count\n"
        ))
        proc = run(good, tmp)
        check("valid sheet passes", proc.returncode == 0, proc.stderr)
        rows = resolved(tmp)
        check("use_for preserved", rows["R"]["use_for"] == "count,stringtie", rows["R"])

        bad = write_sheet(tmp, (
            "sample_id\tfastq_1\tstrandedness\n"
            f"R\t{tmp}/R.fastq.gz\tbackwards\n"
        ), name="bad.tsv")
        proc = run(bad, tmp)
        check("bad strandedness rejected", proc.returncode != 0)


def main() -> int:
    for test in [
        test_explicit_and_inference,
        test_se_with_mate_is_an_error,
        test_pe_without_mate_is_an_error,
        test_mismatched_read_ids_rejected,
        test_duplicate_sample_ids,
        test_missing_file_message,
        test_comments_aliases_and_replicates,
        test_unequal_replicate_lists,
        test_bam_and_fastq_mutually_exclusive,
        test_csv_uses_semicolon_for_replicates,
        test_use_for_and_strandedness_validation,
    ]:
        test()
    print(f"\n{PASS} passed, {FAIL} failed")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
