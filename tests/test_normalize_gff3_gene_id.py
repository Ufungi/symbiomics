#!/usr/bin/env python3
"""Unit tests for bin/normalize_gff3_gene_id.py. No pytest dependency:

    python3 tests/test_normalize_gff3_gene_id.py
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(os.path.dirname(HERE), "bin", "normalize_gff3_gene_id.py")

PASS, FAIL = 0, 0

HIERARCHICAL_GFF3 = """\
##gff-version 3
chr1\tBRAKER\tgene\t1000\t5000\t.\t+\t.\tID=g1
chr1\tBRAKER\tmRNA\t1000\t5000\t.\t+\t.\tID=g1.t1;Parent=g1
chr1\tBRAKER\texon\t1000\t1500\t.\t+\t.\tID=g1.t1.exon1;Parent=g1.t1
chr1\tBRAKER\texon\t3000\t3500\t.\t+\t.\tID=g1.t1.exon2;Parent=g1.t1
chr1\tBRAKER\tCDS\t1000\t1500\t.\t+\t0\tID=g1.t1.CDS1;Parent=g1.t1
chr1\tBRAKER\tgene\t8000\t9000\t.\t-\t.\tID=g2
chr1\tBRAKER\tmRNA\t8000\t9000\t.\t-\t.\tID=g2.t1;Parent=g2
chr1\tBRAKER\texon\t8000\t9000\t.\t-\t.\tID=g2.t1.exon1;Parent=g2.t1
"""

ALREADY_NORMALIZED = """\
##gff-version 3
chr1\tensembl\texon\t100\t200\t.\t+\t.\tgene_id=g1;transcript_id=g1.t1
chr1\tensembl\texon\t300\t400\t.\t+\t.\tgene_id=g1;transcript_id=g1.t1
chr1\tensembl\texon\t500\t600\t.\t+\t.\tgene_id=g2;transcript_id=g2.t1
"""

UNRESOLVABLE = """\
##gff-version 3
chr1\tweird\texon\t100\t200\t.\t+\t.\tID=e1;Parent=missing_parent
chr1\tweird\texon\t300\t400\t.\t+\t.\tID=e2;Parent=also_missing
"""

PINUS_STYLE = """\
##gff-version 3
chr1\tannotation\tgene\t1000\t5000\t.\t+\t.\tID=g1
chr1\tannotation\texon\t1000\t1500\t.\t+\t.\tID=g1
chr1\tannotation\texon\t3000\t3500\t.\t+\t.\tID=g1
chr1\tannotation\tgene\t8000\t9000\t.\t+\t.\tID=g2
chr1\tannotation\texon\t8000\t9000\t.\t+\t.\tID=g2
"""


def check(name: str, condition: bool, detail: str = "") -> None:
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  ok   {name}")
    else:
        FAIL += 1
        print(f"  FAIL {name}  {detail}")


def run_normalize(gff_body: str, tmp: str, *extra: str):
    ann = os.path.join(tmp, "in.gff3")
    out = os.path.join(tmp, "out.gff3")
    report = os.path.join(tmp, "report.json")
    with open(ann, "w") as fh:
        fh.write(gff_body)
    proc = subprocess.run(
        [sys.executable, SCRIPT,
         "--annotation", ann,
         "--out", out,
         "--report", report,
         "--want-feature", "exon",
         *extra],
        capture_output=True, text=True,
    )
    out_body = ""
    if os.path.exists(out):
        with open(out) as fh:
            out_body = fh.read()
    report_body = {}
    if os.path.exists(report):
        with open(report) as fh:
            report_body = json.load(fh)
    return proc, out_body, report_body


def exon_lines(body: str) -> list[str]:
    return [ln for ln in body.splitlines()
            if ln and not ln.startswith("#") and "\texon\t" in ln]


def main() -> int:
    with tempfile.TemporaryDirectory() as tmp1, \
         tempfile.TemporaryDirectory() as tmp2, \
         tempfile.TemporaryDirectory() as tmp3, \
         tempfile.TemporaryDirectory() as tmp4:

        # 1. hierarchical GFF3: Parent chain resolves to gene_id
        proc, out, report = run_normalize(HIERARCHICAL_GFF3, tmp1)
        check("hierarchical exits 0", proc.returncode == 0, proc.stderr)
        exons = exon_lines(out)
        check("3 exon lines", len(exons) == 3, f"got {len(exons)}")
        check("exon1 has gene_id=g1",
              any("gene_id=g1" in ln and "exon1" in ln for ln in exons),
              "\n".join(exons))
        check("exon2 has gene_id=g1",
              any("gene_id=g1" in ln and "exon2" in ln for ln in exons))
        check("g2 exon has gene_id=g2",
              any("gene_id=g2" in ln for ln in exons))
        cds_lines = [ln for ln in out.splitlines() if "\tCDS\t" in ln]
        check("CDS line untouched (not want-feature)",
              len(cds_lines) == 1 and "gene_id=" not in cds_lines[0],
              "\n".join(cds_lines))
        check("report counts 3 covered",
              report.get("n_with_gene_id") == 3, str(report))
        check("report added=3", report.get("n_added") == 3, str(report))

        # 2. already-normalized: passthrough, no rewrite errors
        proc, out, report = run_normalize(ALREADY_NORMALIZED, tmp2)
        check("already-normalized exits 0", proc.returncode == 0, proc.stderr)
        check("already-normalized added=0", report.get("n_added") == 0, str(report))
        check("gene_id preserved",
              out.count("gene_id=g1") == 2 and out.count("gene_id=g2") == 1)

        # 3. unresolvable Parent: coverage 0 -> hard error
        proc, out, report = run_normalize(UNRESOLVABLE, tmp3)
        check("unresolvable exits 1", proc.returncode == 1, proc.stderr)
        check("unresolvable mentions coverage",
              "carry gene_id" in proc.stderr or "coverage" in proc.stderr.lower(),
              proc.stderr)
        check("unresolvable writes no usable output", out == "")

        # 4. Pinus-style: exon ID repeats the gene ID
        proc, out, report = run_normalize(PINUS_STYLE, tmp4)
        check("pinus-style exits 0", proc.returncode == 0, proc.stderr)
        exons = exon_lines(out)
        check("pinus exons get gene_id",
              all("gene_id=g" in ln for ln in exons) and len(exons) == 3,
              "\n".join(exons))

    print(f"\n{PASS} passed, {FAIL} failed")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
