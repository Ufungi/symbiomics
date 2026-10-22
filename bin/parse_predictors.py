#!/usr/bin/env python3
"""Merge SignalP6 / TargetP2 / DeepLoc2 / DeepTMHMM / Phobius / EffectorP /
MEROPS output into one per-protein localisation table.

Every input is optional; a column group is written only for tools that ran.

Secretome call (classical, signal-peptide route):
    n_sp_tools = number of available tools (SignalP6, TargetP2 SP, Phobius,
                 DeepTMHMM SP, DeepLoc2 'Extracellular'/signal) predicting a
                 signal peptide
    secreted   = n_sp_tools >= min(--secretome-min-tools, n tools available)
                 AND no transmembrane helix outside the signal peptide
                 (Phobius TM == 0 and DeepTMHMM class not TM/SP+TM, when run)
The consensus is deliberately a vote, not a single caller: each predictor has
different false-positive behaviour, and agreement is what the fungal
secretome literature (Predector, Sperschneider & Dodds) relies on.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys


def read_ids(path):
    with open(path) as fh:
        return [l[1:].split()[0] for l in fh if l.startswith(">")]


def table(path, sep="\t", comment="#"):
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(comment):
                continue
            yield line.rstrip("\n").split(sep)


def signalp6(p):
    return {f[0].split()[0]: f[1] != "OTHER" for f in table(p) if len(f) > 1}


def targetp2(p):
    out = {}
    for f in table(p):
        if len(f) > 1:
            out[f[0].split()[0]] = f[1]
    return out


def deeploc2(p):
    out = {}
    with open(p, newline="") as fh:
        for r in csv.DictReader(fh):
            out[r["Protein_ID"].split()[0]] = (r.get("Localizations", ""), r.get("Signals", ""))
    return out


def deeptmhmm(p):
    """3line: >id | CLASS / sequence / topology. CLASS in GLOB, SP, TM, SP+TM, BETA."""
    out, last, k = {}, None, 0
    with open(p) as fh:
        for line in (l.rstrip("\n") for l in fh):
            if line.startswith(">"):
                hdr = line[1:].split("|")
                last = hdr[0].strip().split()[0]
                out[last] = [hdr[1].strip() if len(hdr) > 1 else "", 0]
                k = 0
            elif last is not None:
                k += 1
                if k == 2:  # third line of the record: topology string
                    out[last][1] = len(re.findall(r"M+|B+", line))
    return out


def phobius(p):
    out = {}
    for f in table(p, sep=None, comment="SEQ"):
        parts = f if isinstance(f, list) else f.split()
        if len(parts) >= 4:
            out[parts[0]] = (int(parts[1]), parts[2] == "Y")
    return out


def effectorp(p):
    out = {}
    for f in table(p):
        if f:
            out[f[0].split()[0]] = f[-1]
    return out


def merops(p):
    out = {}
    for f in table(p):
        if len(f) >= 2 and f[0] not in out:
            m = re.search(r"\b([ACGMNPSTUI]\d{1,2}[A-Z]?)\b", f[-1] if len(f) > 5 else f[1])
            out[f[0]] = m.group(1) if m else f[1]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--proteome", required=True)
    for k in ("signalp6", "targetp2", "deeploc2", "deeptmhmm", "phobius", "effectorp", "merops"):
        ap.add_argument(f"--{k}")
    ap.add_argument("--secretome-min-tools", type=int, default=2)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    ids = read_ids(a.proteome)
    sp = signalp6(a.signalp6) if a.signalp6 else None
    tp = targetp2(a.targetp2) if a.targetp2 else None
    dl = deeploc2(a.deeploc2) if a.deeploc2 else None
    tm = deeptmhmm(a.deeptmhmm) if a.deeptmhmm else None
    ph = phobius(a.phobius) if a.phobius else None
    ef = effectorp(a.effectorp) if a.effectorp else None
    me = merops(a.merops) if a.merops else None
    sp_sources = [x for x in (sp, tp, ph, tm, dl) if x is not None]
    need = min(a.secretome_min_tools, len(sp_sources))

    cols = ["protein_id"]
    for name, d in (("signalp6", sp), ("targetp2", tp), ("deeploc2", dl), ("deeptmhmm", tm), ("phobius", ph)):
        if d is not None:
            cols.append(name)
    cols += ["n_sp_tools", "n_tm", "secreted"]
    if ef is not None:
        cols.append("effectorp")
    if me is not None:
        cols.append("merops")
    n_sec = 0
    with open(a.out, "w") as o:
        o.write("\t".join(cols) + "\n")
        for pid in ids:
            votes, row, ntm = 0, [pid], 0
            if sp is not None:
                v = sp.get(pid, False); votes += v; row.append("SP" if v else "-")
            if tp is not None:
                v = tp.get(pid, "noTP"); votes += v == "SP"; row.append(v)
            if dl is not None:
                loc, sig = dl.get(pid, ("", ""))
                votes += "Signal peptide" in sig; row.append(loc or "-")
            if tm is not None:
                cls, n = tm.get(pid, ["", 0]); votes += cls.startswith("SP"); ntm = max(ntm, n)
                row.append(f"{cls}:{n}")
            if ph is not None:
                n, s = ph.get(pid, (0, False)); votes += s; ntm = max(ntm, n)
                row.append(f"SP={'Y' if s else 'N'};TM={n}")
            secreted = bool(sp_sources) and votes >= max(need, 1) and (ntm == 0)
            n_sec += secreted
            row += [str(votes), str(ntm), ("yes" if secreted else "no") if sp_sources else "NA"]
            if ef is not None:
                row.append(ef.get(pid, "-") if secreted else "-")
            if me is not None:
                row.append(me.get(pid, "-"))
            o.write("\t".join(row) + "\n")
    sys.stderr.write(f"INFO  ~ [predictors] {len(ids):,} proteins, {n_sec:,} secreted "
                     f"(>= {max(need,1)} of {len(sp_sources)} SP tools, no TM)\n")


if __name__ == "__main__":
    main()
