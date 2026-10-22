"""Offline quantification-first HTML report using only the Python standard library.

All labels are escaped and all report links stay inside the run directory.
Reporting can run before an R engine exists and again after analyses finish.
"""
from __future__ import annotations

import csv
import hashlib
import html
import json
import math
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote

METRICS = {
    "assignment_rate_pct": ("Gene assignment rate", "%", "assigned_count / counting_processed_count × 100; declared physical read/fragment unit"),
    "assignment_counter_rate_pct": ("Native counter assignment (%)", "%", "Assigned / processed native counter events × 100 under declared compatible exclusive semantics; read/fragment unit may be unverified. Counting-stage fraction; sequencing mapping rate is a separate metric."),
    "genes_detected_pct": ("Detected genes / quantified genes", "%", "genes_detected / quantified_gene_count × 100; full gene universe only"),
    "assigned_count": ("Assigned units", "count", "Units assigned by the counting tool; distinct from a gene count sum"),
    "gene_count_sum": ("Gene count sum", "count", "Sum of supplied gene counts; excludes special counters"),
    "genes_detected": ("Detected genes", "genes", "Supplied gene rows with raw count > 0"),
    "quantified_gene_count": ("Quantified genes", "genes", "Number of supplied gene rows before analysis filtering"),
    "mapping_rate_pct": ("Alignment / mapping rate", "%", "Mapped / processed units from a mapping log; distinct from gene assignment"),
    "supplied_gene_detection_pct": ("Detected genes / supplied rows", "%", "Detection within supplied rows; not a whole annotation detection rate"),
}
PALETTE = ["#0f766e", "#2563eb", "#d97706", "#9333ea", "#dc2626", "#64748b", "#be185d", "#4d7c0f"]


def _esc(value):
    return html.escape(str(value if value is not None else ""), quote=True)


def _number(value):
    if value is None or isinstance(value, bool) or str(value).strip().lower() in {"", "na", "nan", "none", "null"}:
        return None
    try:
        result = float(value)
    except (ValueError, TypeError):
        return None
    return result if math.isfinite(result) else None


def _fmt(value):
    value = _number(value)
    if value is None:
        return "NA"
    if value == 0:
        return "0"
    if abs(value) < 0.001 or abs(value) >= 1e9:
        return f"{value:.4e}"
    if value.is_integer():
        return f"{value:,.0f}"
    return f"{value:,.6g}"


def _label(row, column):
    value = row.get(column, "")
    return str(value) if value is not None and str(value).strip() else "(missing)"


def _scope(row, metric=None):
    mapping = metric == "mapping_rate_pct"
    unit = _label(row, "mapping_unit" if mapping else "count_unit")
    denominator = (row.get("mapping_scope") if mapping else
                   row.get("assignment_denominator_scope") or row.get("assignment_scope") or
                   row.get("denominator_scope") or row.get("counting_scope"))
    scope = str(denominator or "unspecified counting scope")
    scope += "; tool=" + _label(row, "quant_tool")
    scope += "; QC=" + str(row.get("qc_level") or "library")
    scope += "; gene universe=" + _label(row, "gene_universe")
    if metric == "assignment_counter_rate_pct":
        scope += "; native counter semantics=" + _label(row, "assignment_mode")
    return _label(row, "assay_id"), unit, scope


def _quantile(values, fraction):
    ordered = sorted(values)
    index = (len(ordered) - 1) * fraction
    low = int(index)
    high = min(low + 1, len(ordered) - 1)
    return ordered[low] + (ordered[high] - ordered[low]) * (index - low)


def _write_tsv(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows({k: "NA" if row.get(k) is None else row.get(k, "") for k in fields} for row in rows)


def _href(run_dir, path):
    """Only existing, local output artifacts can become report links."""
    try:
        target = Path(path)
        if not target.is_absolute():
            target = run_dir / target
        relative = target.resolve().relative_to(run_dir.resolve())
        if not target.exists() or target.is_symlink():
            return None
        return quote(relative.as_posix(), safe="/")
    except (ValueError, OSError, TypeError):
        return None


def _link(run_dir, path, label=None):
    href = _href(run_dir, path)
    return f'<a href="{_esc(href)}">{_esc(label or path)}</a>' if href else _esc(label or path or "")


def _table(rows, fields, limit=None):
    shown = rows if limit is None else rows[:limit]
    if not shown:
        return '<p class="muted">No records.</p>'
    return '<div class="table-scroll"><table><thead><tr>' + ''.join(f'<th>{_esc(k)}</th>' for k in fields) + '</tr></thead><tbody>' + ''.join('<tr>' + ''.join(f'<td>{_esc(_fmt(row.get(k)) if k in METRICS or k in {"counting_processed_count", "n_total", "n_valid", "n_missing", "n_zero", "median", "q1", "q3", "min", "max"} else row.get(k, ""))}</td>' for k in fields) + '</tr>' for row in shown) + '</tbody></table></div>'


def _summary(rows, metric, x, facet, color):
    groups = defaultdict(list)
    for row in rows:
        key = _scope(row, metric) + ((_label(row, facet) if facet else "all"), _label(row, x), (_label(row, color) if color else "all"))
        groups[key].append(row)
    result = []
    for key, members in sorted(groups.items()):
        values = [n for row in members if (n := _number(row.get(metric))) is not None]
        assay, unit, scope, facet_value, x_value, color_value = key
        result.append({
            "metric": metric, "assay_id": assay, "count_unit": unit, "denominator_scope": scope,
            "group_column": x, "group_value": x_value, "facet_column": facet, "facet_value": facet_value,
            "color_column": color, "color_value": color_value, "n_total": len(members), "n_valid": len(values),
            "n_missing": len(members) - len(values), "n_zero": sum(v == 0 for v in values),
            "n_biological_units": len({str(r["biological_unit_id"]) for r in members if r.get("biological_unit_id")}),
            "median": _quantile(values, .5) if values else None,
            "q1": _quantile(values, .25) if values else None, "q3": _quantile(values, .75) if values else None,
            "min": min(values) if values else None, "max": max(values) if values else None,
        })
    return result


def _svg_scatter(rows, metric, x, color, scale, title, order=""):
    """Grouped sample dots with medians/IQR. On log axes zero has a separate lane."""
    values = [(r, _number(r.get(metric))) for r in rows]
    valid = [(r, n) for r, n in values if n is not None and n >= 0]
    groups = list(dict.fromkeys(_label(r, x) for r in rows))
    ordered = [part for part in str(order).split("|") if part in groups]
    groups = ordered + [g for g in groups if g not in ordered]
    colors = sorted({_label(r, color) if color else "samples" for r in rows})
    color_map = {value: PALETTE[i % len(PALETTE)] for i, value in enumerate(colors)}
    width = max(760, min(5000, 180 + 92 * len(groups)))
    left, right, top, bottom = 90, width - 35, 72, 327
    legend_lines = math.ceil(max(len(colors), 1) / 3)
    height = 490 + 24 * legend_lines
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-label="{_esc(title)}">',
           '<style>text{font:13px system-ui,sans-serif;fill:#334155}.grid{stroke:#e2e8f0;stroke-width:1}.axis{stroke:#64748b}</style>',
           f'<rect width="{width}" height="{height}" fill="white"/>',
           f'<text x="22" y="28" style="font-size:18px;font-weight:600">{_esc(title)}</text>']
    positives = [n for _, n in valid if n > 0]
    is_log = scale == "log10"
    if is_log:
        low = math.floor(math.log10(min(positives))) - 1 if positives else -1
        high = math.ceil(math.log10(max(positives))) if positives else 0
        if high <= low:
            high = low + 1
        ticks = [low + (high - low) * i / min(6, high - low) for i in range(min(6, high - low) + 1)]
        coord = lambda n: bottom + 30 if n == 0 else bottom - (math.log10(n) - low) / (high - low) * (bottom - top)
        svg.append(f'<line x1="{left}" x2="{right}" y1="{bottom + 30}" y2="{bottom + 30}" stroke="#f59e0b" stroke-dasharray="4 4"/>')
        svg.append(f'<text x="{left-10}" y="{bottom+34}" text-anchor="end">true 0</text>')
        tick_values = [(coord(10 ** t), _fmt(10 ** t)) for t in ticks]
    else:
        high = max((n for _, n in valid), default=1) or 1
        high *= 1.08
        coord = lambda n: bottom - n / high * (bottom - top)
        tick_values = [(coord(high * i / 5), _fmt(high * i / 5)) for i in range(6)]
    for y, label in tick_values:
        svg.append(f'<line x1="{left}" x2="{right}" y1="{y:.2f}" y2="{y:.2f}" class="grid"/>')
        svg.append(f'<text x="{left-10}" y="{y+4:.2f}" text-anchor="end">{_esc(label)}</text>')
    unit = METRICS.get(metric, (metric, "count", ""))[1]
    svg.append(f'<text x="18" y="{top-19}">{_esc(unit)}{(" · log10" if is_log else " · linear")}</text>')
    spacing = (right - left) / max(1, len(groups))
    group_x = {g: left + (i + .5) * spacing for i, g in enumerate(groups)}
    for group in groups:
        gx = group_x[group]
        members = [(r, n) for r, n in valid if _label(r, x) == group]
        # Box summaries are on the raw values, transformed only for display.
        for ci, clabel in enumerate(colors):
            nums = [n for r, n in members if (_label(r, color) if color else "samples") == clabel]
            if not nums:
                continue
            bx = gx + (ci - (len(colors) - 1) / 2) * min(14, spacing / max(len(colors) + 1, 2))
            q1, median, q3 = (_quantile(nums, f) for f in (.25, .5, .75))
            svg.extend([f'<line x1="{bx:.2f}" x2="{bx:.2f}" y1="{coord(q1):.2f}" y2="{coord(q3):.2f}" stroke="{color_map[clabel]}" stroke-width="7" opacity=".35"/>',
                        f'<line x1="{bx-6:.2f}" x2="{bx+6:.2f}" y1="{coord(median):.2f}" y2="{coord(median):.2f}" stroke="{color_map[clabel]}" stroke-width="2"/>'])
        for row, value in members:
            label = _label(row, "sample_id")
            seed = int(hashlib.sha256(label.encode()).hexdigest()[:8], 16) / 0xffffffff - .5
            px = gx + seed * min(58, spacing * .65)
            py = coord(value)
            clabel = _label(row, color) if color else "samples"
            tooltip = _esc(f"{label} | {group} | {clabel}: {value:.12g} {unit}")
            if value == 0:
                svg.append(f'<path d="M {px-3:.2f},{py-3:.2f} l 6,6 m -6,0 l 6,-6" stroke="{color_map[clabel]}" stroke-width="1.8"><title>{tooltip} (true zero)</title></path>')
            else:
                svg.append(f'<circle cx="{px:.2f}" cy="{py:.2f}" r="3.2" fill="{color_map[clabel]}" opacity=".75"><title>{tooltip}</title></circle>')
        label_y = bottom + (62 if is_log else 32)
        svg.append(f'<text x="{gx:.2f}" y="{label_y}" text-anchor="end" transform="rotate(-25 {gx:.2f} {label_y})">{_esc(group)}</text>')
        svg.append(f'<text x="{gx:.2f}" y="{label_y+47}" text-anchor="middle" style="font-size:11px">n={len(members)}/{sum(_label(r,x)==group for r in rows)}</text>')
    legend_y = height - 34 - 24 * (legend_lines - 1)
    for i, label in enumerate(colors):
        px = 24 + (i % 3) * (width - 40) / 3
        py = legend_y + (i // 3) * 24
        svg.append(f'<circle cx="{px}" cy="{py-4}" r="4" fill="{color_map[label]}"/><text x="{px+12}" y="{py}">{_esc(label)}</text>')
    svg.append('</svg>')
    return '\n'.join(svg)


def _read_tsv(path):
    try:
        with path.open(encoding="utf-8-sig", newline="") as stream:
            return list(csv.DictReader(stream, delimiter="\t"))
    except (OSError, UnicodeError):
        return []


def _svg_coordinates(rows, xcol, ycol, color, color_levels, axis_labels, title, bounds):
    """Two-dimensional scores; all facets share axes and metadata colors."""
    width, height = 1000, max(650, 130 + 24 * len(color_levels))
    left, right, top, bottom = 90, 750, 70, 550
    xmin, xmax, ymin, ymax = bounds
    xp = lambda n: left + (n - xmin) / (xmax - xmin) * (right - left)
    yp = lambda n: bottom - (n - ymin) / (ymax - ymin) * (bottom - top)
    color_map = {label: PALETTE[i] if i < len(PALETTE) else f"hsl({i * 137.507 % 360:.2f},58%,40%)" for i, label in enumerate(color_levels)}
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-label="{_esc(title)}">',
           '<style>text{font:13px system-ui,sans-serif;fill:#334155}.grid{stroke:#e2e8f0;stroke-width:1}</style>',
           f'<rect width="{width}" height="{height}" fill="white"/>',
           f'<text x="24" y="30" style="font-size:18px;font-weight:600">{_esc(title)}</text>']
    for index in range(6):
        xv, yv = xmin + (xmax - xmin) * index / 5, ymin + (ymax - ymin) * index / 5
        px, py = xp(xv), yp(yv)
        svg.extend([f'<line x1="{px:.2f}" x2="{px:.2f}" y1="{top}" y2="{bottom}" class="grid"/>',
                    f'<text x="{px:.2f}" y="{bottom+24}" text-anchor="middle">{_esc(_fmt(xv))}</text>',
                    f'<line x1="{left}" x2="{right}" y1="{py:.2f}" y2="{py:.2f}" class="grid"/>',
                    f'<text x="{left-10}" y="{py+4:.2f}" text-anchor="end">{_esc(_fmt(yv))}</text>'])
    if xmin < 0 < xmax:
        svg.append(f'<line x1="{xp(0):.2f}" x2="{xp(0):.2f}" y1="{top}" y2="{bottom}" stroke="#94a3b8" stroke-dasharray="4 4"/>')
    if ymin < 0 < ymax:
        svg.append(f'<line x1="{left}" x2="{right}" y1="{yp(0):.2f}" y2="{yp(0):.2f}" stroke="#94a3b8" stroke-dasharray="4 4"/>')
    for row in rows:
        label = _label(row, color) if color else "samples"
        tooltip = _esc(f"{row['sample_id']} | {label} | {xcol}={row[xcol]}, {ycol}={row[ycol]}")
        svg.append(f'<circle cx="{xp(float(row[xcol])):.3f}" cy="{yp(float(row[ycol])):.3f}" r="5" fill="{color_map[label]}" stroke="white" opacity=".85"><title>{tooltip}</title></circle>')
    for index, label in enumerate(color_levels):
        py = 83 + index * 24
        svg.append(f'<circle cx="780" cy="{py-4}" r="5" fill="{color_map[label]}"/><text x="792" y="{py}">{_esc(label)}</text>')
    svg.extend([f'<text x="{(left+right)/2}" y="{bottom+60}" text-anchor="middle">{_esc(axis_labels[0])}</text>',
                f'<text x="23" y="{(top+bottom)/2}" transform="rotate(-90 23 {(top+bottom)/2})" text-anchor="middle">{_esc(axis_labels[1])}</text>',
                f'<text x="{left}" y="{height-16}" style="font-size:12px">n={len(rows)} selected observations; linear axes; shared limits across facets; hover for sample ID</text>', '</svg>'])
    return '\n'.join(svg)


def _analysis_plots(run_dir, analyses, requests):
    """Draw requested projections from engine coordinates without recalculating PCA/MDS."""
    from urllib.parse import unquote
    known = {str(a.get("analysis_id", "")) for a in analyses}
    outcomes = []
    for index, request in enumerate(requests):
        if request.get("scope") != "analysis":
            continue
        aid, metric = str(request.get("analysis_id", "")), request.get("metric", "")
        outcome = {"plot_id": request.get("plot_id", str(index)), "analysis_id": aid, "metric": metric, "status": "blocked", "reason": "", "path": "", "_request_index": request.get("_request_index", index)}
        try:
            if aid not in known:
                raise ValueError(f"unknown analysis_id {aid!r}")
            if metric not in {"pca", "mds"}:
                raise ValueError(f"unsupported analysis metric {metric!r}; use pca or mds")
            if request.get("scale") not in {None, "", "linear"}:
                raise ValueError("PCA/MDS coordinate plots require scale=linear")
            root = run_dir / "analyses" / aid
            if not root.resolve().is_relative_to((run_dir / "analyses").resolve()):
                raise ValueError("analysis path is outside this run")
            scores_path = root / "exploratory" / f"{metric}_scores.tsv"
            samples_path = root / "input" / "samples.tsv"
            if not scores_path.exists():
                outcome["reason"] = f"Coordinates unavailable: {scores_path.relative_to(run_dir)}; analysis may be pending, validation-only or unable to normalize."
                outcomes.append(outcome)
                continue
            if not samples_path.exists():
                raise ValueError("analysis input/samples.tsv is missing; metadata cannot be joined")
            scores, metadata = _read_tsv(scores_path), _read_tsv(samples_path)
            if not scores or not metadata:
                raise ValueError("coordinate or sample metadata table is empty/unreadable")
            ids, mids = [r.get("sample_id") for r in scores], [r.get("sample_id") for r in metadata]
            if None in ids or None in mids or "" in ids or "" in mids or len(set(ids)) != len(ids) or len(set(mids)) != len(mids):
                raise ValueError("coordinate and metadata tables need unique nonempty sample_id values")
            lookup = {r["sample_id"]: r for r in metadata}
            if set(ids) - set(lookup):
                raise ValueError("coordinate sample_id values are missing from recorded metadata")
            color, facet = request.get("color") or "", request.get("facet") or ""
            for column in (color, facet):
                if column and not all(column in row for row in metadata):
                    raise ValueError(f"missing metadata column {column!r}")
            xcol, ycol = ("PC1", "PC2") if metric == "pca" else ("MDS1", "MDS2")
            joined = []
            for score in scores:
                if _number(score.get(xcol)) is None or _number(score.get(ycol)) is None:
                    raise ValueError(f"nonfinite or missing {xcol}/{ycol} coordinate for {score['sample_id']}")
                # Coordinate axes come exclusively from engine output; metadata joins by ID.
                joined.append({**lookup[score["sample_id"]], xcol: score[xcol], ycol: score[ycol]})
            color_levels = sorted({_label(row, color) if color else "samples" for row in joined})
            order = [unquote(value, encoding="utf-8", errors="strict") for value in str(request.get("order") or "").split("|") if value]
            ordered = list(dict.fromkeys(level for level in order if level in color_levels))
            color_levels = ordered + [level for level in color_levels if level not in ordered]
            xs, ys = [float(r[xcol]) for r in joined], [float(r[ycol]) for r in joined]
            # Equal coordinate units per pixel, shared by all facet panels.
            midpoint_x, midpoint_y = (min(xs)+max(xs))/2, (min(ys)+max(ys))/2
            units_per_pixel = max((max(xs)-min(xs))/660, (max(ys)-min(ys))/480, 1e-12) * 1.15
            bounds = (midpoint_x-units_per_pixel*330, midpoint_x+units_per_pixel*330,
                      midpoint_y-units_per_pixel*240, midpoint_y+units_per_pixel*240)
            axis_labels = [xcol, ycol]
            if metric == "pca":
                fractions = {r.get("component"): _number(r.get("variance_fraction")) for r in _read_tsv(root / "exploratory" / "pca_variance.tsv")}
                axis_labels = [f"PC{i} ({fractions[str(i)]*100:.1f}% variance)" if fractions.get(str(i)) is not None else f"PC{i}" for i in (1, 2)]
            panels = defaultdict(list)
            for row in joined:
                panels[_label(row, facet) if facet else "all"].append(row)
            dest = root / "exploratory" / "configured_plots"
            dest.mkdir(parents=True, exist_ok=True)
            settings = {"plot_id": outcome["plot_id"], "metric": metric, "source": scores_path.relative_to(run_dir).as_posix(), "metadata_source": samples_path.relative_to(run_dir).as_posix(), "sample_join": "sample_id", "axes": [xcol, ycol], "x_setting_ignored": request.get("x") or "", "color": color, "facet": facet, "color_order": color_levels, "scale": "linear", "bounds": bounds, "sample_count": len(joined), "panels": []}
            for panel_index, (facet_value, members) in enumerate(sorted(panels.items())):
                stem = f"plot_{index:02d}_{metric}_{panel_index:02d}"
                title = f"{aid}: {metric.upper()}" + (f" — {facet}={facet_value}" if facet else "")
                image = dest / (stem + ".svg")
                image.write_text(_svg_coordinates(members, xcol, ycol, color, color_levels, axis_labels, title, bounds), encoding="utf-8")
                data = dest / (stem + ".tsv")
                _write_tsv(data, members, list(dict.fromkeys(["sample_id", xcol, ycol] + [col for col in (color, facet, "biological_unit_id", "metadata_status", "biological_unit_status") if col])))
                settings["panels"].append({"facet": facet_value, "n": len(members), "image": image.relative_to(run_dir).as_posix(), "data": data.relative_to(run_dir).as_posix()})
            (dest / f"plot_{index:02d}_{metric}_settings.json").write_text(json.dumps(settings, ensure_ascii=False, indent=2), encoding="utf-8")
            outcome.update(status="ok", reason=f"Recorded {xcol}/{ycol}; color={color or 'none'}, facet={facet or 'none'}. x setting is ignored; coordinates remain unchanged.", path=dest.relative_to(run_dir).as_posix())
        except (ValueError, OSError, UnicodeError) as exc:
            outcome.update(status="failed", reason=str(exc))
        outcomes.append(outcome)
    _write_tsv(run_dir / "analysis_plot_status.tsv", outcomes, ["plot_id", "analysis_id", "metric", "status", "reason", "path"])
    return outcomes


def _validate_plot_requests(plots, rows):
    """Keep one outcome per requested table row, including invalid and duplicate IDs."""
    from collections import Counter
    ids = Counter(str(p.get("plot_id") or "") for p in plots)
    outcomes, accepted = [], []
    for index, request in enumerate(plots):
        plot_id = str(request.get("plot_id") or "")
        scope, metric = request.get("scope", ""), request.get("metric", "")
        outcome = {"request_row": index + 1, "plot_id": plot_id, "scope": scope,
                   "analysis_id": request.get("analysis_id") or "", "metric": metric,
                   "status": "blocked", "reason": "Plot has not been rendered yet", "path": ""}
        errors = []
        if not plot_id:
            errors.append("plot_id is required")
        if ids[plot_id] > 1:
            errors.append(f"duplicate plot_id {plot_id!r}; all duplicated requests are rejected")
        if scope not in {"quant_qc", "analysis"}:
            errors.append(f"unsupported plot scope {scope!r}; use quant_qc or analysis")
        elif scope == "quant_qc":
            if metric not in METRICS:
                errors.append(f"unsupported QC metric {metric!r}")
            scale = request.get("scale") or "linear"
            if scale not in {"linear", "log10"}:
                errors.append(f"unsupported QC scale {scale!r}; use linear or log10")
            columns = [request.get("x") or "study_id", request.get("color") or "", request.get("facet") or ""]
            missing = [col for col in columns if col and rows and not any(col in row for row in rows)]
            if missing:
                errors.append("missing metadata columns: " + ", ".join(missing))
        if errors:
            outcome.update(status="failed", reason="; ".join(errors))
        else:
            accepted.append({**request, "_request_index": index})
        outcomes.append(outcome)
    return accepted, outcomes


def render_report(run_dir: Path, metrics: list[dict], analyses: list[dict], statuses: list[dict], provenance: list | dict, plots: list[dict] | None = None) -> Path:
    """Write report.html, grouped descriptive summaries and SVG/data sidecars.

    Metadata may be flattened in metrics or nested in a ``metadata`` mapping.
    Every requested plot has an outcome in plot_status.tsv; invalid requests fail.
    The caller owns sample_metrics.tsv, metric_definitions.tsv, unassigned_long.tsv.
    """
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    qc = run_dir / "00_quant_qc"
    plot_dir = qc / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    rows = [{**(r.get("metadata") if isinstance(r.get("metadata"), dict) else {}), **r} for r in metrics]
    defaults = [
        {"plot_id": "assignment_rate", "scope": "quant_qc", "metric": "assignment_rate_pct", "x": "study_id", "scale": "log10"},
        {"plot_id": "detection_rate", "scope": "quant_qc", "metric": "genes_detected_pct", "x": "study_id", "scale": "linear"},
        {"plot_id": "gene_count_sum", "scope": "quant_qc", "metric": "gene_count_sum", "x": "study_id", "scale": "log10"},
        {"plot_id": "assigned_count", "scope": "quant_qc", "metric": "assigned_count", "x": "study_id", "scale": "log10"},
    ]
    if any(_number(r.get("assignment_counter_rate_pct")) is not None for r in rows):
        defaults.insert(1, {"plot_id": "native_counter_assignment", "scope": "quant_qc", "metric": "assignment_counter_rate_pct", "x": "study_id", "scale": "log10"})
    if any(r.get("gene_universe") != "full" for r in rows):
        detection_index = next(i for i, item in enumerate(defaults) if item["metric"] == "genes_detected_pct") + 1
        defaults.insert(detection_index, {"plot_id": "supplied_detection", "scope": "quant_qc", "metric": "supplied_gene_detection_pct", "x": "study_id", "scale": "linear"})
    if any(_number(r.get("mapping_rate_pct")) is not None for r in rows):
        defaults.append({"plot_id": "mapping_rate", "scope": "quant_qc", "metric": "mapping_rate_pct", "x": "study_id", "scale": "linear"})
    accepted_requests, plot_outcomes = _validate_plot_requests(plots or [], rows)
    requests = defaults + [p for p in accepted_requests if p.get("scope", "quant_qc") == "quant_qc"]
    summaries, figures, unavailable = [], [], []
    for index, request in enumerate(requests):
        metric = request.get("metric", "assignment_rate_pct")
        if metric not in METRICS:
            unavailable.append(f"{request.get('plot_id', index)}: unsupported QC metric {metric}")
            continue
        x, color, facet = request.get("x") or "study_id", request.get("color") or "", request.get("facet") or ""
        scale = request.get("scale") or "linear"
        if scale not in {"linear", "log10"}:
            unavailable.append(f"{request.get('plot_id', index)}: unsupported scale {scale}")
            continue
        missing_cols = [c for c in (x, color, facet) if c and rows and not any(c in r for r in rows)]
        if missing_cols and index >= len(defaults):
            unavailable.append(f"{request.get('plot_id', index)}: missing metadata columns {', '.join(missing_cols)}")
            continue
        summaries.extend(_summary(rows, metric, x, facet, color))
        panels = defaultdict(list)
        rendered_panels, missing_panels = 0, 0
        for row in rows:
            panels[_scope(row, metric) + ((_label(row, facet) if facet else "all"),)].append(row)
        for panel_index, (key, members) in enumerate(sorted(panels.items())):
            if not any(_number(r.get(metric)) is not None for r in members):
                unavailable.append(f"{METRICS[metric][0]} — {key[0]} / {key[1]} / {key[3]}: NA for all {len(members)} samples; see metric reasons.")
                missing_panels += 1
                continue
            rendered_panels += 1
            name = f"qc_{index:02d}_{panel_index:02d}"
            title = METRICS[metric][0]
            if metric in {"supplied_gene_detection_pct", "genes_detected_pct"}:
                denominators = {_number(row.get("quantified_gene_count")) for row in members}
                if None not in denominators and len(denominators) == 1:
                    denominator = next(iter(denominators))
                    universe_label = "supplied" if metric == "supplied_gene_detection_pct" else "quantified"
                    title = f"Detected out of {_fmt(denominator)} {universe_label} genes (%)"
            subtitle = f"assay={key[0]}; count unit={key[1]}; {key[2]}" + (f"; {facet}={key[3]}" if facet else "")
            image_path = plot_dir / f"{name}.svg"
            image_path.write_text(_svg_scatter(members, metric, x, color, scale, title, request.get("order", "")), encoding="utf-8")
            data_path = plot_dir / f"{name}.tsv"
            sidecar_fields = list(dict.fromkeys(["sample_id", "assay_id", "count_unit", "denominator_scope", "assigned_count", "counting_processed_count", "assignment_counter_status", "assignment_counter_reason", metric, x] + [v for v in (color, facet, "biological_unit_id", "status", "reason") if v]))
            _write_tsv(data_path, [{**r, "denominator_scope": _scope(r, metric)[2]} for r in members], sidecar_fields)
            native_note = " Counting-stage native-counter fraction under declared exclusive semantics; physical read/fragment unit may be unverified." if metric == "assignment_counter_rate_pct" else ""
            figures.append(f'<figure><figcaption><strong>{_esc(title)}</strong><br><span class="muted">{_esc(subtitle)}</span></figcaption><div class="plot-scroll"><img src="{_esc(_href(run_dir,image_path))}" alt="{_esc(title + '; ' + subtitle)}" loading="lazy"></div><p class="muted">x: {_esc(x)}; scale: {_esc(scale)}. Dots = samples; × = true zero; thick line = IQR; short line = median. n = valid/all observations.{native_note} {_link(run_dir,data_path,"Plot data TSV")}</p></figure>')
        if "_request_index" in request:
            outcome = plot_outcomes[request["_request_index"]]
            if rendered_panels:
                outcome.update(status="ok", reason=f"Rendered {rendered_panels} QC panels; {missing_panels} panels have no available values.", path="00_quant_qc/plots")
            else:
                outcome.update(status="blocked", reason=f"No available {metric} values for the requested QC plot; inspect metric availability reasons.")
    # Deduplicate defaults and user-requested summaries by their complete semantic key.
    summaries = list({tuple((k, str(v)) for k, v in s.items()): s for s in summaries}.values())
    summary_fields = ["metric", "assay_id", "count_unit", "denominator_scope", "group_column", "group_value", "facet_column", "facet_value", "color_column", "color_value", "n_total", "n_valid", "n_missing", "n_zero", "n_biological_units", "median", "q1", "q3", "min", "max"]
    _write_tsv(qc / "group_summary.tsv", summaries, summary_fields)
    input_rows = [r for r in rows if r.get("qc_level") != "technical_aggregate"]
    assignment_values = [_number(r.get("assignment_rate_pct")) for r in input_rows]
    zero_assigned = sum(_number(r.get("assigned_count")) == 0 for r in input_rows)
    zero_gene = sum(_number(r.get("gene_count_sum")) == 0 for r in input_rows)
    cards = [("Input observations", len(input_rows)), ("Assigned count = 0", zero_assigned), ("Gene count sum = 0", zero_gene), ("Known-unit assignment rate unavailable", sum(v is None for v in assignment_values)), ("Small positive rate (<0.01%)", sum(v is not None and 0 < v < .01 for v in assignment_values))]
    if any("assignment_counter_rate_pct" in r for r in input_rows):
        counter_values = [_number(r.get("assignment_counter_rate_pct")) for r in input_rows]
        cards[-1] = ("Small positive native-counter rate (<0.01%)", sum(v is not None and 0 < v < .01 for v in counter_values))
        cards.insert(4, ("Native-counter rate unavailable", sum(v is None for v in counter_values)))
    card_html = ''.join(f'<div class="card"><b>{n:,}</b><span>{_esc(label)}</span></div>' for label, n in cards)
    scope_rows = []
    scope_groups = defaultdict(list)
    for row in rows:
        scope_groups[_scope(row)].append(row)
    for (assay, unit, scope), members in sorted(scope_groups.items()):
        complete = lambda metric: sum(_number(r.get(metric)) for r in members) if unit not in {"unknown", "(missing)"} and all(_number(r.get(metric)) is not None for r in members) else None
        scope_rows.append({"assay": assay, "unit": unit, "counting scope": scope, "observations": len(members), "assigned total": _fmt(complete("assigned_count")), "processed total": _fmt(complete("counting_processed_count")), "gene count total": _fmt(complete("gene_count_sum"))})
    fields = ["sample_id", "assay_id", "qc_level", "count_unit", "gene_count_sum", "assigned_count", "counting_processed_count", "assignment_rate_pct", "assignment_counter_rate_pct", "assignment_counter_status", "assignment_counter_reason", "genes_detected", "quantified_gene_count", "genes_detected_pct", "mapping_rate_pct", "assignment_status", "assignment_reason", "status", "reason"]
    reason_columns = sorted({k for row in rows for k in row if k.endswith("_reason") or k.endswith("_scope") or k in {"gene_universe", "provenance_level", "assignment_mode", "include_analysis", "exclusion_reason", "metadata_status"}})
    extra = [k for k in reason_columns if k not in fields]
    sample_table = _table(rows, fields + extra)
    status_html = '<table><thead><tr><th>Stage</th><th>Analysis</th><th>Status</th><th>Reason</th><th>Output</th></tr></thead><tbody>' + ''.join(f'<tr><td>{_esc(s.get("stage",""))}</td><td>{_esc(s.get("analysis_id",""))}</td><td><span class="badge {_esc(s.get("status",""))}">{_esc(s.get("status",""))}</span></td><td>{_esc(s.get("reason",""))}</td><td>{_link(run_dir,s.get("path","")) if s.get("path") else ""}</td></tr>' for s in statuses) + '</tbody></table>'
    for outcome in _analysis_plots(run_dir, analyses, accepted_requests):
        request_index = outcome.pop("_request_index")
        plot_outcomes[request_index].update(outcome)
    analysis_plot_statuses = [outcome for outcome in plot_outcomes if outcome["scope"] == "analysis"]
    _write_tsv(run_dir / "analysis_plot_status.tsv", analysis_plot_statuses, ["plot_id", "analysis_id", "metric", "status", "reason", "path"])
    _write_tsv(run_dir / "plot_status.tsv", plot_outcomes, ["request_row", "plot_id", "scope", "analysis_id", "metric", "status", "reason", "path"])
    unavailable.extend(f"{s['plot_id']}: {s['reason']}" for s in plot_outcomes if s["status"] != "ok")
    plot_status_html = _table(plot_outcomes, ["request_row", "plot_id", "scope", "analysis_id", "metric", "status", "reason"]) if plot_outcomes else '<p>No custom plots were requested; default QC and engine figures are shown.</p>' 
    analysis_sections = []
    for analysis in analyses:
        aid = str(analysis.get("analysis_id", ""))
        root = run_dir / "analyses" / aid
        safe = root.resolve().is_relative_to((run_dir / "analyses").resolve())
        artifacts = sorted(root.rglob("*")) if safe and root.exists() else []
        links = ''.join(f'<li>{_link(run_dir,p,p.relative_to(root).as_posix())}</li>' for p in artifacts if p.is_file() and p.suffix.lower() in {".tsv", ".csv", ".json", ".txt", ".pdf", ".png", ".svg", ".rds"})
        previews = ''.join(f'<figure><figcaption>{_esc(p.relative_to(root).as_posix())}</figcaption><img class="analysis-image" src="{_esc(_href(run_dir,p))}" alt="{_esc(p.name)}" loading="lazy"></figure>' for p in artifacts if p.is_file() and p.suffix.lower() in {".png", ".svg"} and _href(run_dir,p))
        mode = analysis.get("analysis_mode", "inferential")
        engine_rows = []
        if safe and (root / "status.json").exists():
            try:
                engine_status = json.loads((root / "status.json").read_text(encoding="utf-8"))
                stages = engine_status.get("stages", {})
                if isinstance(stages, dict):
                    engine_rows = [{"stage": stage, **outcome} for stage, outcome in stages.items() if isinstance(outcome, dict)]
            except (ValueError, OSError):
                engine_rows = [{"stage": "status", "status": "unavailable", "reason": "Could not read engine status.json"}]
        pdf_links = ''.join(f'<li>{_link(run_dir,p,p.relative_to(root).as_posix())}</li>' for p in artifacts if p.is_file() and p.suffix.lower() == ".pdf")
        pdf_section = f'<p>Analysis figures (PDF):</p><ul>{pdf_links}</ul>' if pdf_links else ""
        configured_status = [r for r in analysis_plot_statuses if r["analysis_id"] == aid]
        configured_table = _table(configured_status, ["plot_id", "metric", "status", "reason"]) if configured_status else ""
        pending = '<p class="muted">No analysis artifacts are available yet; inspect stage outcomes below.</p>' if not artifacts else ""
        analysis_sections.append(f'<article><h3>{_esc(aid)} <small>{_esc(mode)}</small></h3>{_table([analysis],list(analysis))}{_table(engine_rows,["stage","status","reason"]) if engine_rows else ""}{configured_table}{pending}{previews}{pdf_section}<details><summary>All analysis artifacts</summary><ul>{links}</ul></details></article>')
    definitions = [{"metric": k, "unit": v[1], "definition": v[2]} for k, v in METRICS.items()]
    unassigned = _read_tsv(qc / "unassigned_long.tsv")
    unassigned_panels = defaultdict(list)
    metric_lookup = {(r.get("sample_id"), r.get("qc_level", "library")): r for r in rows}
    for item in unassigned:
        metric_row = metric_lookup.get((item.get("sample_id"), item.get("qc_level", "library")), {})
        datum = {**metric_row, **item, "unassigned_count": item.get("count")}
        unassigned_panels[_scope(metric_row)].append(datum)
    for index, (key, members) in enumerate(sorted(unassigned_panels.items())):
        path = plot_dir / f"unassigned_{index:02d}.svg"
        data_path = path.with_suffix(".tsv")
        title = "Original counting categories — sample counts"
        path.write_text(_svg_scatter(members, "unassigned_count", "category", "", "log10", title), encoding="utf-8")
        _write_tsv(data_path, members, ["sample_id", "assay_id", "tool", "category", "count", "unit", "composition_valid", "qc_level"])
        figures.append(f'<figure><figcaption><strong>{_esc(title)}</strong><br><span class="muted">{_esc("; ".join(key))}</span></figcaption><div class="plot-scroll"><img src="{_esc(_href(run_dir,path))}" alt="Original counting category counts" loading="lazy"></div><p class="muted">Each dot is a sample/category count. Values use a log10 axis; true zero is shown separately. Categories retain original tool definitions and may include Assigned. {_link(run_dir,data_path,"Plot data TSV")}</p></figure>')
    unassigned_html = '<p>Original counting categories are preserved. Missing categories are not assumed to be zero. Overlapping categories must not be interpreted as a 100% composition.</p>' + _table(unassigned, list(unassigned[0]) if unassigned else [])
    qc_links = ' · '.join(_link(run_dir,qc/name,name) for name in ["sample_metrics.tsv", "metric_definitions.tsv", "group_summary.tsv", "unassigned_long.tsv"] if (qc/name).exists())
    warnings = '<ul>' + ''.join(f'<li>{_esc(w)}</li>' for w in unavailable) + '</ul>' if unavailable else ''
    provenance_text = json.dumps(provenance, ensure_ascii=False, indent=2, default=str)
    failed = [s for s in statuses if s.get("status") in {"failed", "blocked"}]
    run_state = "Partial / failed stages — inspect status records" if failed else ("Completed" if isinstance(provenance, dict) and provenance.get("status") == "ok" else "In progress — available stage outcomes below")
    doc = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<meta http-equiv="Content-Security-Policy" content="default-src 'none'; img-src 'self' data:; style-src 'unsafe-inline'; base-uri 'none'; form-action 'none'">
<title>Symbiomics — {_esc(run_dir.name)}</title><style>
:root{{font-family:system-ui,-apple-system,sans-serif;color:#172b3a;background:#f1f5f9;line-height:1.55}}body{{margin:0}}main{{max-width:1220px;margin:auto;padding:32px 24px 70px}}header{{border-bottom:3px solid #0f766e;margin-bottom:28px}}h1{{font-size:32px;margin:0}}h2{{font-size:24px;margin-top:36px}}h3{{font-size:20px}}small,.muted{{color:#64748b}}a{{color:#0f766e}}nav{{display:flex;gap:22px;flex-wrap:wrap;margin:18px 0}}.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:12px}}.card,figure,article,.notice{{background:white;border:1px solid #dbe4ea;border-radius:10px;padding:18px}}.card b{{display:block;font-size:29px;color:#0f766e}}.card span{{font-size:13px}}.notice{{border-left:4px solid #d97706}}figure{{margin:18px 0}}figcaption{{margin-bottom:10px}}.plot-scroll,.table-scroll{{overflow:auto;max-width:100%}}.plot-scroll img{{display:block;width:100%;min-width:760px}}.analysis-image{{max-width:100%;height:auto}}table{{border-collapse:collapse;font-size:12px;width:100%;background:white}}th,td{{padding:9px 11px;border-bottom:1px solid #e2e8f0;text-align:left;vertical-align:top;white-space:nowrap}}td:last-child{{white-space:normal;min-width:100px}}th{{background:#e8f0f2;color:#23404a}}details{{margin:16px 0}}summary{{cursor:pointer;font-weight:600}}pre{{white-space:pre-wrap;overflow-wrap:anywhere;background:#e8eef3;padding:16px;border-radius:8px;font-size:12px}}.badge{{padding:2px 8px;border-radius:20px;background:#e2e8f0}}.ok{{background:#ccfbf1;color:#115e59}}.failed,.blocked{{background:#fee2e2;color:#991b1b}}article{{margin:18px 0}}footer{{margin-top:36px;color:#64748b;font-size:12px}}@media print{{body{{background:white}}details{{display:block}}.plot-scroll img{{min-width:0}}main{{padding:0}}}}
</style></head><body><main><header><p class="muted">SYMBIOMICS · BULK RNA-SEQ POST-ANALYSIS</p><h1>Quantification comes first</h1><p>Run {_esc(run_dir.name)} · {_esc(run_state)}</p><nav><a href="#qc">Quantification QC</a><a href="#samples">Samples &amp; reasons</a><a href="#analyses">Analyses</a><a href="#status">Stage status</a><a href="#provenance">Provenance</a></nav></header>
<section id="qc"><h2>All input observations</h2><div class="cards">{card_html}</div><p>{qc_links}</p><p class="notice">NA means unavailable, not zero. Small positive values retain scientific notation and log-scale positions; true zeros use × markers and a separate lane on log plots. Counts and rates are summarized separately by assay, count unit, counting scope, tool and QC level. Cards count input observations; technical aggregates are displayed separately. n counts observations; biological units are recorded separately in group_summary.tsv.</p>
{_table(scope_rows,["assay","unit","counting scope","observations","assigned total","processed total","gene count total"])}<p class="muted">Totals are NA when any constituent numerator or denominator is unavailable or the physical count unit is unknown. Native-counter proportions remain sample-specific; unknown physical units are not pooled into count totals. Technical partitions may appear as multiple input observations; consult technical merge mappings before counting biological replicates.</p>
{''.join(figures)}<details><summary>Unavailable metrics / plot requests ({len(unavailable)})</summary>{warnings}</details><details><summary>Metric definitions</summary>{_table(definitions,["metric","unit","definition"])}</details><details><summary>Grouped descriptive statistics</summary>{_table(summaries,summary_fields)}</details><details><summary>Unassigned categories</summary>{unassigned_html}</details></section>
<section id="samples"><h2>Sample metrics and availability</h2><p>Low signal and zero observations remain in this table. Metadata declarations and per-metric reasons are retained; report labels do not certify unverified provenance.</p>{sample_table}</section>
<section id="analyses"><h2>Analysis designs and results</h2><p class="notice">Exploratory PCA describes the selected samples. DEG requires a confirmed experimental design and independent biological units. Within-target relative expression and the fraction of reads assigned to the target answer different questions; neither quantity alone establishes organism abundance, colonization or causation.</p>{''.join(analysis_sections) if analysis_sections else '<p>No statistical analysis was configured for this run. Quantification QC is available above.</p>'}</section>
<section id="status"><h2>Requested plot outcomes</h2><p>{_link(run_dir,'plot_status.tsv','All plot request outcomes TSV')}</p>{plot_status_html}<h2>Stage outcomes</h2><div class="table-scroll">{status_html}</div><p>Failed or blocked stages do not erase available QC or completed analysis outputs. Skipped stages state why they were not run.</p></section>
<section id="provenance"><h2>Inputs and execution provenance</h2><p>{_link(run_dir,'manifest.json','Manifest JSON')} · {_link(run_dir,'stage_status.tsv','Stage status TSV')}</p><details><summary>Recorded provenance</summary><pre>{_esc(provenance_text)}</pre></details></section>
<footer>Standalone local report. No JavaScript, fonts, analytics or other network resources are loaded. Numeric summaries are descriptive; no automatic group significance tests are added.</footer></main></body></html>'''
    path = run_dir / "report.html"
    path.write_text(doc, encoding="utf-8")
    return path
