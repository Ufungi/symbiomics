#!/usr/bin/env python3
"""Turn a DIAMOND blastp-vs-Swiss-Prot outfmt6 table into per-query best-hit
descriptions, by parsing the standard UniProt FASTA header directly.

Why not full UPIMAPI here: UPIMAPI's rich column set (GO/EC/Pathway/KEGG/Pfam)
comes from its ID-mapping step, which needs a downloaded ID-mapping resource
this lab does not have on disk yet (Phase B, see docs/functional_annotation.md).
Every UniProt FASTA entry already carries its description in the header
(`>sp|ACCESSION|NAME DESCRIPTION OS=... OX=... GN=... PE=... SV=...`), so a
DIAMOND best hit plus a header parse gets the single most load-bearing field --
the product description -- with zero extra downloads. This is Phase A's
Swiss-Prot leg; full UPIMAPI with GO/EC/Pathway columns is Phase B.
"""
from __future__ import annotations

import argparse
import re
import sys

HEADER_RE = re.compile(
    r"^(?P<db>sp|tr)\|(?P<acc>[^|]+)\|(?P<entry>\S+)\s+(?P<desc>.*?)"
    r"(?:\s+OS=(?P<os>.*?))?(?:\s+OX=\d+)?(?:\s+GN=\S+)?(?:\s+PE=\d+)?(?:\s+SV=\d+)?$"
)


def index_fasta_headers(path: str) -> dict[str, dict]:
    index = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            header = line[1:].rstrip("\n")
            token = header.split(None, 1)[0]
            acc = token.split("|")[1] if "|" in token else token
            match = HEADER_RE.match(header)
            if match:
                index[acc] = {
                    "entry": match.group("entry"),
                    "description": match.group("desc").strip(),
                    "organism": (match.group("os") or "").strip(),
                }
            else:
                index[acc] = {"entry": token, "description": header, "organism": ""}
    return index


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--diamond-tsv", required=True,
                    help="DIAMOND outfmt6: qseqid sseqid pident length evalue bitscore")
    ap.add_argument("--swissprot-fasta", required=True)
    ap.add_argument("--min-pident", type=float, default=30.0)
    ap.add_argument("--max-evalue", type=float, default=1e-5)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    headers = index_fasta_headers(args.swissprot_fasta)
    sys.stderr.write(f"INFO  ~ [swissprot] indexed {len(headers):,} Swiss-Prot headers\n")

    best: dict[str, tuple] = {}
    with open(args.diamond_tsv) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            qid, sid, pident, _length, evalue, bitscore = fields[:6]
            pident, evalue, bitscore = float(pident), float(evalue), float(bitscore)
            if pident < args.min_pident or evalue > args.max_evalue:
                continue
            if qid not in best or bitscore > best[qid][1]:
                best[qid] = (sid, bitscore, pident, evalue)

    n_hit = 0
    with open(args.out, "w") as out:
        out.write("query_id\tswissprot_accession\tswissprot_entry\tpident\tevalue\tbitscore\tdescription\torganism\n")
        for qid, (sid, bitscore, pident, evalue) in sorted(best.items()):
            acc = sid.split("|")[1] if "|" in sid else sid
            info = headers.get(acc, {"entry": sid, "description": "", "organism": ""})
            out.write(
                f"{qid}\t{acc}\t{info['entry']}\t{pident:.1f}\t{evalue:.2e}\t{bitscore:.1f}\t"
                f"{info['description']}\t{info['organism']}\n"
            )
            n_hit += 1

    sys.stderr.write(f"INFO  ~ [swissprot] {n_hit:,} queries with a hit >= {args.min_pident}% id\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
