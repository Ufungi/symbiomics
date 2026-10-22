"""Multi-organism read composition from HTSeq special counters.

When reads are aligned separately to several genomes (host + symbiont/pathogen),
each library has one HTSeq raw output per genome. The first (primary) organism's
counters define the input read-pair total (assigned + all `__*` categories);
every organism contributes aligned pairs (total - `__not_aligned`) and assigned
pairs. Optionally one organism defines a signal group (reads per million input
pairs against a threshold) and a QC pass flag on primary aligned pairs.

Signal groups describe read abundance only; they are not an infection rate or
dose, and no downstream statistic should treat them as one.
"""
import csv
from pathlib import Path

try:
    from .importers import write_tsv
except ImportError:
    from importers import write_tsv


def read_htseq_counters(path):
    assigned = special = not_aligned = 0
    with Path(path).open(encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            gene, _, n = line.rstrip("\n").partition("\t")
            n = int(n.split("\t")[0])
            if gene.startswith("__"):
                special += n
                if gene == "__not_aligned":
                    not_aligned = n
            else:
                assigned += n
    return dict(assigned=assigned, total=assigned + special, not_aligned=not_aligned)


def read_libraries(path, base):
    rows = list(csv.DictReader(Path(path).open(encoding="utf-8", newline=""), delimiter="\t"))
    need = {"sample_id", "organism", "file"}
    if not rows or not need <= set(rows[0]):
        raise ValueError(f"{path} needs columns {sorted(need)}")
    for r in rows:
        f = Path(r["file"])
        r["file"] = f if f.is_absolute() else Path(base) / f
        if not r["file"].exists():
            raise ValueError(f"missing HTSeq file for {r['sample_id']}/{r['organism']}: {r['file']}")
    return rows


def compute(rows, primary, signal=None, signal_ppm=1000.0, min_aligned_pairs=0.0):
    by = {}
    for r in rows:
        by.setdefault(r["sample_id"], {})[r["organism"]] = read_htseq_counters(r["file"])
    organisms = list(dict.fromkeys(r["organism"] for r in rows))
    if primary not in organisms:
        raise ValueError(f"primary organism '{primary}' not in libraries table ({organisms})")
    if signal and signal not in organisms:
        raise ValueError(f"signal organism '{signal}' not in libraries table ({organisms})")
    out = []
    for sid, d in by.items():
        if primary not in d:
            raise ValueError(f"{sid}: no {primary} counters; the primary organism defines the input total")
        total = d[primary]["total"]
        row = dict(sample_id=sid, input_pairs=total)
        for org in organisms:
            c = d.get(org)
            aligned = (c["total"] - c["not_aligned"]) if c else None
            row[f"{org}_aligned_pairs"] = aligned
            row[f"{org}_aligned_pct"] = round(100 * aligned / total, 4) if c and total else None
            row[f"{org}_assigned_pairs"] = c["assigned"] if c else None
            row[f"{org}_assigned_ppm"] = round(1e6 * c["assigned"] / total, 4) if c and total else None
        row["qc_pass"] = bool(row[f"{primary}_aligned_pairs"] >= min_aligned_pairs)
        if signal:
            ppm = row[f"{signal}_assigned_ppm"]
            row["signal_group"] = "NA" if ppm is None else ("high" if ppm >= signal_ppm else "low")
        out.append(row)
    return organisms, out


COLORS = ["#5B8C5A", "#C97B2D", "#4C78A8", "#B279A2", "#E45756"]


def stacked_svg(organisms, rows, path):
    """Stacked horizontal bars: % of input pairs aligned per organism; remainder grey."""
    h, left, width = 16, 150, 520
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{left + width + 190}" height="{40 + h * len(rows) + 30 * len(organisms)}" font-family="sans-serif" font-size="11">']
    for i, r in enumerate(rows):
        y = 20 + i * h
        parts.append(f'<text x="{left - 6}" y="{y + 11}" text-anchor="end">{r["sample_id"]}</text>')
        x = 0.0
        for k, org in enumerate(organisms):
            pct = r[f"{org}_aligned_pct"] or 0
            w = min(pct, 100 - x) / 100 * width
            parts.append(f'<rect x="{left + x / 100 * width:.1f}" y="{y}" width="{w:.1f}" height="{h - 3}" fill="{COLORS[k % len(COLORS)]}"/>')
            x += min(pct, 100 - x)
        parts.append(f'<rect x="{left + x / 100 * width:.1f}" y="{y}" width="{(100 - x) / 100 * width:.1f}" height="{h - 3}" fill="#d9d9d9"/>')
    ly = 30 + len(rows) * h
    for k, org in enumerate(organisms + ["unaligned / other"]):
        col = COLORS[k % len(COLORS)] if k < len(organisms) else "#d9d9d9"
        parts.append(f'<rect x="{left}" y="{ly + k * 14}" width="10" height="10" fill="{col}"/><text x="{left + 16}" y="{ly + k * 14 + 9}">{org}</text>')
    parts.append("</svg>")
    Path(path).write_text("\n".join(parts), encoding="utf-8")


def run(libraries, out_dir, primary, signal=None, signal_ppm=1000.0, min_aligned_pairs=0.0):
    out = Path(out_dir)
    rows = read_libraries(libraries, Path(libraries).resolve().parent)
    organisms, table = compute(rows, primary, signal, signal_ppm, min_aligned_pairs)
    write_tsv(out / "composition.tsv", table)
    stacked_svg(organisms, table, out / "composition.svg")
    return out / "composition.tsv"
