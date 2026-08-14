#!/usr/bin/env python3
"""Assign one product description per protein, from whichever functional
sources ran, in a fixed priority order. No AHRD: this is a direct best-hit
ladder rather than a weighted-blast synthesis step.

Priority, first non-empty wins:
    1. funannotate2 product           (fungi branch only, Phase B / F4)
    2. Swiss-Prot best hit            (DIAMOND, this pipeline's Phase A)
    3. eggNOG-mapper Description / Preferred_name
    4. InterProScan best signature description  (Phase B)
    5. "hypothetical protein"

Every source is optional -- pass only the ones that actually ran. A protein
with no source at all still gets a row, so the denominator for %hypothetical
is always every input protein, not just the ones something matched.
"""
from __future__ import annotations

import argparse
import csv
import sys

PLACEHOLDER_DESCRIPTIONS = {
    "", "-", "n/a", "na", "unknown", "unnamed protein product",
    "uncharacterized protein", "hypothetical protein",
}


def read_ids(path: str) -> list[str]:
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def read_swissprot(path: str) -> dict[str, dict]:
    out = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            desc = row.get("description", "").strip()
            if desc and desc.lower() not in PLACEHOLDER_DESCRIPTIONS:
                out[row["query_id"]] = {
                    "product": desc,
                    "accession": row.get("swissprot_accession", ""),
                    "pident": row.get("pident", ""),
                }
    return out


def read_funannotate2(path: str) -> dict[str, dict]:
    out = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            desc = (row.get("product") or row.get("Product") or "").strip()
            gid = row.get("gene_id") or row.get("ID") or row.get("query_id")
            if gid and desc and desc.lower() not in PLACEHOLDER_DESCRIPTIONS:
                out[gid] = {"product": desc}
    return out


def read_eggnog(path: str) -> dict[str, dict]:
    """emapper.py .annotations: '##'-prefixed metadata, one '#query' header row."""
    out = {}
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if header is None or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            row = dict(zip(header, fields))
            qid = row.get("query", "")
            if not qid:
                continue
            desc = (row.get("Description") or "").strip()
            pref = (row.get("Preferred_name") or "").strip()
            product = desc if desc and desc != "-" else pref
            entry = {
                "go": row.get("GOs", ""), "ec": row.get("EC", ""),
                "kegg_ko": row.get("KEGG_ko", ""), "cog": row.get("COG_category", ""),
            }
            if product and product.lower() not in PLACEHOLDER_DESCRIPTIONS and product != "-":
                entry["product"] = product
            out[qid] = entry
    return out


def read_interpro(path: str) -> dict[str, dict]:
    """IPS TSV: protein_accession, ..., signature_description, ..."""
    out = {}
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            qid, desc = fields[0], fields[5].strip()
            if desc and desc.lower() not in PLACEHOLDER_DESCRIPTIONS and desc != "-":
                if qid not in out:
                    out[qid] = {"product": desc}
    return out


def read_dbcan(path: str) -> dict[str, dict]:
    """dbCAN3 overview.txt: Gene ID, HMMER, eCAMI, DIAMOND, #ofTools -- the
    column set narrows if RUN_DBCAN was invoked with fewer than all three
    callers. dbCAN never names a plain-English product, so it doesn't enter
    the priority ladder; it only contributes a CAZy family call and counts as
    one more source of evidence for the consensus tag."""
    out = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        gid_col = idx.get("Gene ID", 0)
        tools_col = idx.get("#ofTools")
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= gid_col:
                continue
            gid = fields[gid_col]
            families = set()
            for caller in ("HMMER", "eCAMI", "DIAMOND"):
                i = idx.get(caller)
                if i is None or i >= len(fields):
                    continue
                val = fields[i].strip()
                if val and val != "-":
                    families.add(val.split("+")[0].split("(")[0])
            n_tools = fields[tools_col].strip() if tools_col is not None and tools_col < len(fields) else str(len(families))
            if families or (n_tools not in ("", "0")):
                out[gid] = {"family": ";".join(sorted(families)), "n_tools": n_tools}
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--proteome", required=True, help="FASTA -- defines the full ID set")
    ap.add_argument("--funannotate2", default=None)
    ap.add_argument("--swissprot", default=None)
    ap.add_argument("--eggnog", default=None)
    ap.add_argument("--interpro", default=None)
    ap.add_argument("--dbcan", default=None, help="dbCAN overview.txt (Phase A)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    all_ids = read_ids(args.proteome)
    if not all_ids:
        sys.stderr.write(f"ERROR ~ [products] no sequences found in {args.proteome}\n")
        return 1

    f2 = read_funannotate2(args.funannotate2) if args.funannotate2 else {}
    sp = read_swissprot(args.swissprot) if args.swissprot else {}
    eg = read_eggnog(args.eggnog) if args.eggnog else {}
    ip = read_interpro(args.interpro) if args.interpro else {}
    dc = read_dbcan(args.dbcan) if args.dbcan else {}

    source_counts = {"funannotate2": 0, "swissprot": 0, "eggnog": 0, "interpro": 0, "none": 0}

    with open(args.out, "w") as out:
        out.write(
            "protein_id\tproduct\tsource\tgo_terms\tec_number\tkegg_ko\t"
            "swissprot_hit\tswissprot_pident\tdbcan_family\tdbcan_tools\t"
            "sources_with_evidence\tn_sources\n"
        )
        for pid in all_ids:
            go = eg.get(pid, {}).get("go", "")
            ec = eg.get(pid, {}).get("ec", "")
            ko = eg.get(pid, {}).get("kegg_ko", "")
            sp_hit = sp.get(pid, {}).get("accession", "")
            sp_pid = sp.get(pid, {}).get("pident", "")
            dbcan_family = dc.get(pid, {}).get("family", "")
            dbcan_tools = dc.get(pid, {}).get("n_tools", "")

            if pid in f2:
                product, source = f2[pid]["product"], "funannotate2"
            elif pid in sp:
                product, source = sp[pid]["product"], "swissprot"
            elif "product" in eg.get(pid, {}):
                product, source = eg[pid]["product"], "eggnog"
            elif pid in ip:
                product, source = ip[pid]["product"], "interpro"
            else:
                product, source = "hypothetical protein", "none"

            # Evidence tags are independent of the ladder above: a protein can
            # carry evidence from several sources at once (e.g. swissprot +
            # eggnog + dbcan all hit the same gene) even though only one of
            # them supplies the product name. This is the per-gene "different
            # programs agree" view Blast2GO shows as tags.
            evidence = [name for name, hit in (
                ("funannotate2", pid in f2), ("swissprot", pid in sp),
                ("eggnog", pid in eg), ("dbcan", pid in dc), ("interpro", pid in ip),
            ) if hit]

            source_counts[source] += 1
            out.write(
                f"{pid}\t{product}\t{source}\t{go}\t{ec}\t{ko}\t{sp_hit}\t{sp_pid}\t"
                f"{dbcan_family}\t{dbcan_tools}\t{','.join(evidence)}\t{len(evidence)}\n"
            )

    total = len(all_ids)
    sys.stderr.write(f"INFO  ~ [products] {total:,} proteins:\n")
    for source, count in source_counts.items():
        if count:
            sys.stderr.write(f"        {source:14s} {count:6,} ({count/total*100:5.1f}%)\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
