#!/usr/bin/env python3
"""Render assign_products.py's products.tsv as a self-contained HTML report.

The free, in-pipeline alternative to Blast2GO's annotation table: one static
HTML file, sortable/filterable client-side, no server and no external CDN
(this lab's compute nodes are offline) -- just the TSV data embedded as JSON
plus vanilla JS. The per-gene "tags" column is the same idea as Blast2GO's
consensus tags and dbCAN's own #ofTools column: which independent programs
found evidence for this protein, at a glance.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys

SOURCE_COLORS = {
    "funannotate2": "#7c5cbf",
    "swissprot": "#2b7a78",
    "eggnog": "#c0562f",
    "dbcan": "#2f6fc0",
    "interpro": "#a3358a",
}


def load_rows(path: str) -> list[dict]:
    rows = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            tags = [t for t in row.get("sources_with_evidence", "").split(",") if t]
            rows.append({
                "id": row.get("protein_id", ""),
                "product": row.get("product", ""),
                "source": row.get("source", ""),
                "tags": tags,
                "go": row.get("go_terms", ""),
                "ec": row.get("ec_number", ""),
                "ko": row.get("kegg_ko", ""),
                "dbcanFamily": row.get("dbcan_family", ""),
                "dbcanTools": row.get("dbcan_tools", ""),
                "spHit": row.get("swissprot_hit", ""),
                "spPident": row.get("swissprot_pident", ""),
            })
    return rows


def summarize(rows: list[dict]) -> dict:
    total = len(rows)
    by_source = {}
    tag_counts = {}
    for row in rows:
        by_source[row["source"]] = by_source.get(row["source"], 0) + 1
        for tag in row["tags"]:
            tag_counts[tag] = tag_counts.get(tag, 0) + 1
    hypothetical = by_source.get("none", 0)
    return {
        "total": total,
        "hypotheticalFrac": round(hypothetical / total, 4) if total else 0.0,
        "bySource": by_source,
        "tagCounts": tag_counts,
    }


HTML_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Functional annotation report</title>
<style>
:root {
    --bg: #ffffff; --panel: #f6f7f9; --border: #dde1e6; --text: #1b1f24;
    --muted: #5b6470; --accent: #2b6fb0; --row-alt: #fafbfc;
}
@media (prefers-color-scheme: dark) {
    :root:not([data-theme="light"]) {
        --bg: #14161a; --panel: #1c1f24; --border: #30343b; --text: #e7e9ec;
        --muted: #9099a6; --accent: #6ea8e0; --row-alt: #191c21;
    }
}
* { box-sizing: border-box; }
body {
    margin: 0; background: var(--bg); color: var(--text);
    font: 14px/1.5 -apple-system, "Segoe UI", Helvetica, Arial, sans-serif;
}
header {
    padding: 1.25rem 1.5rem; border-bottom: 1px solid var(--border);
    display: flex; flex-wrap: wrap; gap: 1.5rem; align-items: baseline;
}
header h1 { font-size: 1.1rem; margin: 0; }
.stat { font-size: 0.85rem; color: var(--muted); }
.stat b { color: var(--text); font-variant-numeric: tabular-nums; }
.tagcount { display: inline-flex; align-items: center; gap: 0.3rem; margin-right: 0.75rem; }
.dot { width: 0.6em; height: 0.6em; border-radius: 50%; display: inline-block; }
.controls {
    padding: 0.75rem 1.5rem; display: flex; gap: 0.75rem; align-items: center;
    border-bottom: 1px solid var(--border); background: var(--panel);
}
#filter {
    flex: 1; max-width: 28rem; padding: 0.4rem 0.6rem; border: 1px solid var(--border);
    border-radius: 6px; background: var(--bg); color: var(--text); font-size: 0.9rem;
}
.pagebar { display: flex; gap: 0.5rem; align-items: center; margin-left: auto; color: var(--muted); font-size: 0.85rem; }
.pagebar button {
    border: 1px solid var(--border); background: var(--bg); color: var(--text);
    border-radius: 5px; padding: 0.25rem 0.6rem; cursor: pointer;
}
.pagebar button:disabled { opacity: 0.4; cursor: default; }
main { overflow-x: auto; padding: 0 1.5rem 2rem; }
table { border-collapse: collapse; width: 100%; min-width: 900px; }
th, td { padding: 0.4rem 0.6rem; border-bottom: 1px solid var(--border); text-align: left; vertical-align: top; }
th { position: sticky; top: 0; background: var(--panel); cursor: pointer; user-select: none; white-space: nowrap; }
th:after { content: ""; opacity: 0.4; margin-left: 0.3rem; }
th[data-dir="asc"]:after { content: "▲"; opacity: 1; }
th[data-dir="desc"]:after { content: "▼"; opacity: 1; }
tbody tr:nth-child(even) { background: var(--row-alt); }
.mono { font-family: ui-monospace, "SFMono-Regular", Menlo, Consolas, monospace; font-size: 0.85em; }
.tags { display: flex; flex-wrap: wrap; gap: 0.25rem; }
.tag {
    color: #fff; font-size: 0.72rem; padding: 0.08rem 0.45rem; border-radius: 999px;
    white-space: nowrap;
}
.hypothetical { color: var(--muted); font-style: italic; }
</style>
</head>
<body>
<header>
    <h1>Functional annotation report</h1>
    <span class="stat"><b id="statTotal"></b> proteins</span>
    <span class="stat"><b id="statHypo"></b> hypothetical</span>
    <span class="stat" id="statTags"></span>
</header>
<div class="controls">
    <input id="filter" type="text" placeholder="Filter by protein id, product, or tag...">
    <div class="pagebar">
        <button id="prevPage">&larr; prev</button>
        <span id="pageLabel"></span>
        <button id="nextPage">next &rarr;</button>
    </div>
</div>
<main>
<table>
<thead><tr>
    <th data-key="id">protein_id</th>
    <th data-key="product">product</th>
    <th data-key="source">source</th>
    <th data-key="tags">evidence</th>
    <th data-key="dbcanFamily">CAZy family</th>
    <th data-key="go">GO</th>
    <th data-key="ec">EC</th>
    <th data-key="ko">KEGG KO</th>
    <th data-key="spPident">Swiss-Prot %id</th>
</tr></thead>
<tbody id="rows"></tbody>
</table>
</main>
<script>
const ROWS = __ROWS_JSON__;
const SUMMARY = __SUMMARY_JSON__;
const SOURCE_COLORS = __COLORS_JSON__;
const PAGE_SIZE = 100;

document.getElementById('statTotal').textContent = SUMMARY.total.toLocaleString();
document.getElementById('statHypo').textContent =
    (SUMMARY.bySource.none || 0).toLocaleString() + ' (' + (SUMMARY.hypotheticalFrac * 100).toFixed(1) + '%)';
document.getElementById('statTags').innerHTML = Object.entries(SUMMARY.tagCounts)
    .sort((a, b) => b[1] - a[1])
    .map(([tag, n]) => `<span class="tagcount"><span class="dot" style="background:${SOURCE_COLORS[tag] || '#888'}"></span>${tag} <b>${n.toLocaleString()}</b></span>`)
    .join('');

let sortKey = null, sortDir = 1, page = 0;
let filtered = ROWS;

function applyFilter() {
    const q = document.getElementById('filter').value.trim().toLowerCase();
    filtered = !q ? ROWS : ROWS.filter(r =>
        r.id.toLowerCase().includes(q) ||
        r.product.toLowerCase().includes(q) ||
        r.tags.some(t => t.toLowerCase().includes(q))
    );
    page = 0;
    render();
}

function applySort(key) {
    if (sortKey === key) { sortDir = -sortDir; } else { sortKey = key; sortDir = 1; }
    document.querySelectorAll('th[data-key]').forEach(th => {
        th.removeAttribute('data-dir');
        if (th.dataset.key === sortKey) th.setAttribute('data-dir', sortDir === 1 ? 'asc' : 'desc');
    });
    filtered = filtered.slice().sort((a, b) => {
        const av = Array.isArray(a[key]) ? a[key].length : (a[key] || '');
        const bv = Array.isArray(b[key]) ? b[key].length : (b[key] || '');
        return av < bv ? -sortDir : av > bv ? sortDir : 0;
    });
    page = 0;
    render();
}

function esc(s) {
    return String(s).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
}

function render() {
    const start = page * PAGE_SIZE;
    const slice = filtered.slice(start, start + PAGE_SIZE);
    const tbody = document.getElementById('rows');
    tbody.innerHTML = slice.map(r => `
        <tr>
            <td class="mono">${esc(r.id)}</td>
            <td class="${r.source === 'none' ? 'hypothetical' : ''}">${esc(r.product)}</td>
            <td>${esc(r.source)}</td>
            <td><div class="tags">${r.tags.map(t => `<span class="tag" style="background:${SOURCE_COLORS[t] || '#888'}">${esc(t)}</span>`).join('')}</div></td>
            <td class="mono">${esc(r.dbcanFamily)}</td>
            <td class="mono">${esc(r.go)}</td>
            <td class="mono">${esc(r.ec)}</td>
            <td class="mono">${esc(r.ko)}</td>
            <td>${esc(r.spPident)}</td>
        </tr>`).join('');
    const pages = Math.max(1, Math.ceil(filtered.length / PAGE_SIZE));
    document.getElementById('pageLabel').textContent =
        `${filtered.length.toLocaleString()} rows · page ${page + 1} / ${pages}`;
    document.getElementById('prevPage').disabled = page === 0;
    document.getElementById('nextPage').disabled = page >= pages - 1;
}

document.getElementById('filter').addEventListener('input', applyFilter);
document.getElementById('prevPage').addEventListener('click', () => { if (page > 0) { page--; render(); } });
document.getElementById('nextPage').addEventListener('click', () => {
    if ((page + 1) * PAGE_SIZE < filtered.length) { page++; render(); }
});
document.querySelectorAll('th[data-key]').forEach(th => th.addEventListener('click', () => applySort(th.dataset.key)));

render();
</script>
</body>
</html>
"""


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--products", required=True, help="assign_products.py output TSV")
    ap.add_argument("--out", default="annotation_report.html")
    args = ap.parse_args()

    rows = load_rows(args.products)
    if not rows:
        sys.stderr.write(f"ERROR ~ [report] {args.products} has no data rows\n")
        return 1

    summary = summarize(rows)
    html = (
        HTML_TEMPLATE
        .replace("__ROWS_JSON__", json.dumps(rows, separators=(",", ":")))
        .replace("__SUMMARY_JSON__", json.dumps(summary))
        .replace("__COLORS_JSON__", json.dumps(SOURCE_COLORS))
    )
    with open(args.out, "w") as fh:
        fh.write(html)

    sys.stderr.write(f"INFO  ~ [report] wrote {args.out} ({len(rows):,} proteins)\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
