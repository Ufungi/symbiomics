#!/usr/bin/env python3
"""Merge per-sample TrimGalore trimming reports into a single TSV."""
import re, csv, sys, os, glob

rows = []
for rpt in sorted(glob.glob('*.trimming_report.txt')):
    sample = rpt.split('.trimming_report')[0]
    txt = open(rpt).read()

    def extract(pattern, text, cast=str):
        m = re.search(pattern, text)
        return cast(m.group(1).replace(',', '')) if m else ''

    total      = extract(r'Total reads processed:\s+(\S+)', txt, int)
    adapters   = extract(r'Reads with adapters:\s+(\S+)', txt, int)
    passed     = extract(r'Reads written \(passing filters\):\s+(\S+)', txt, int)
    total_bp   = extract(r'Total basepairs processed:\s+(\S+)', txt, int)
    quality    = extract(r'Quality-trimmed:\s+(\S+)', txt, int)
    written_bp = extract(r'Total written \(filtered\):\s+(\S+)', txt, int)
    rows.append([sample, total, adapters, passed, total_bp, quality, written_bp])

with open('trim_summary.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow([
        'sample_id', 'total_reads', 'reads_with_adapters',
        'reads_passing_filters', 'total_bp', 'quality_trimmed_bp',
        'total_written_bp'
    ])
    for r in sorted(rows, key=lambda x: x[0]):
        w.writerow(r)
