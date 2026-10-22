"""Count adapters and auditable quantification QC; no third-party dependencies.

import_inputs(inputs, samples, base_dir, out_dir) returns assays, metrics,
unassigned, provenance, warnings, errors. Assay compatibility errors preserve QC
and leave counts_path empty. Structural input errors raise ValueError.
write_tsv(path, rows, fieldnames=None) writes None as NA without rounding numbers.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import defaultdict
from decimal import Decimal, InvalidOperation
from pathlib import Path

HTSEQ_COUNTERS = ("no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique")
FORMATS = {"htseq_matrix", "featurecounts_matrix", "htseq_raw", "featurecounts_raw", "salmon_quant", "generic_raw_matrix"}


def write_tsv(path, rows, fieldnames=None):
    """Export dictionaries, preserving integer precision and missing-vs-zero."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    fields = list(fieldnames or dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: "NA" if row.get(k) is None else row.get(k, "") for k in fields})
    return str(path.resolve())


def number(value, context="count", integer=False):
    try:
        value = Decimal(str(value))
    except (InvalidOperation, ValueError):
        raise ValueError(f"{context}: invalid numeric value {value!r}") from None
    if not value.is_finite() or value < 0:
        raise ValueError(f"{context}: expected finite nonnegative count")
    if value == value.to_integral_value():
        return int(value)
    if integer:
        raise ValueError(f"{context}: fractional value in integer count mode")
    result = float(value)
    if not math.isfinite(result) or (result == 0 and value != 0):
        raise ValueError(f"{context}: count outside supported numeric range")
    return result


def total(values):
    values = list(values)
    return sum(values) if all(isinstance(v, int) for v in values) else math.fsum(values)


def percentage(numerator, denominator):
    return 100 * (numerator / denominator) if numerator is not None and denominator and denominator > 0 else None


def resolve(base, value):
    path = Path(value).expanduser()
    return path if path.is_absolute() else Path(base) / path


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def read_table(path, comments=False):
    path = Path(path)
    with path.open(encoding="utf-8-sig", newline="") as handle:
        lines = (line for line in handle if not (comments and line.startswith("#")))
        reader = csv.reader(lines, delimiter="," if path.suffix.lower() == ".csv" else "\t")
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError(f"{path}: empty table") from None
        if not header or any(not h for h in header) or len(set(header)) != len(header):
            raise ValueError(f"{path}: missing or duplicate header")
        rows = []
        for line, values in enumerate(reader, 2):
            if not values or not any(values):
                continue
            if len(values) != len(header):
                raise ValueError(f"{path}:{line}: expected {len(header)} fields, got {len(values)}")
            rows.append(dict(zip(header, values)))
    return header, rows


def _matrix(path, fmt, integer):
    header, rows = read_table(path, comments=fmt == "featurecounts_raw")
    if fmt == "featurecounts_raw":
        if header[:6] != ["Geneid", "Chr", "Start", "End", "Strand", "Length"]:
            raise ValueError(f"{path}: invalid featureCounts header")
        columns = header[6:]
    else:
        if header[0] not in ("gene_id", "Geneid"):
            raise ValueError(f"{path}: first column must be gene_id")
        columns = header[1:]
    if not columns or not rows:
        raise ValueError(f"{path}: empty count matrix")
    values = {column: {} for column in columns}
    for row in rows:
        gene = row[header[0]]
        if not gene or gene in values[columns[0]]:
            raise ValueError(f"{path}: empty or duplicate gene ID {gene!r}")
        if gene.startswith("__"):
            raise ValueError(f"{path}: HTSeq special counter {gene!r} must not be in a gene matrix")
        for column in columns:
            values[column][gene] = number(row[column], f"{path}: {gene}/{column}", integer)
    return values


def _htseq(path, integer):
    genes, special = {}, {}
    with Path(path).open(encoding="utf-8-sig") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 2:
                raise ValueError(f"{path}:{line_no}: expected two HTSeq columns")
            key, value = fields
            target = special if key.startswith("__") else genes
            key = key[2:] if key.startswith("__") else key
            if not key or key in target:
                raise ValueError(f"{path}:{line_no}: duplicate/empty feature {key!r}")
            target[key] = number(value, f"{path}:{line_no}", integer if target is genes else True)
    if not genes:
        raise ValueError(f"{path}: no gene counts")
    return genes, special


def _salmon(path, mapping_path):
    _, rows = read_table(path)
    if not rows or not {"Name", "Length", "EffectiveLength", "TPM", "NumReads"}.issubset(rows[0]):
        raise ValueError(f"{path}: invalid/empty Salmon quant.sf")
    if not mapping_path:
        raise ValueError(f"{path}: Salmon requires explicit transcript_id/gene_id mapping")
    header, mappings = read_table(mapping_path)
    if not {"transcript_id", "gene_id"}.issubset(header):
        raise ValueError(f"{mapping_path}: expected transcript_id and gene_id headers")
    txmap = {}
    for row in mappings:
        tid, gid = row["transcript_id"], row["gene_id"]
        if not tid or not gid or tid in txmap:
            raise ValueError(f"{mapping_path}: duplicate or empty mapping {tid!r}")
        txmap[tid] = gid
    seen, genes = set(), defaultdict(list)
    for row in rows:
        tid = row["Name"]
        if tid in seen or tid not in txmap:
            raise ValueError(f"{path}: duplicate/unmapped transcript {tid!r}")
        seen.add(tid)
        for field in ("Length", "EffectiveLength", "TPM"):
            number(row[field], f"{path}: {tid}/{field}")
        genes[txmap[tid]].append(number(row["NumReads"], f"{path}: {tid}/NumReads"))
    return {gene: total(values) for gene, values in genes.items()}


def _legacy_row(path, sample_key):
    header, rows = read_table(path)
    keys = [k for k in ("sample_id_original", "sample", "sample_id") if k in header]
    # Piden carries intentional display-label swaps. Match only the canonical
    # source-ID column; searching aliases together would match two samples.
    if not keys:
        raise ValueError(f"{path}: no supported statistics sample-ID column")
    found = [row for row in rows if row[keys[0]] == sample_key]
    if len(found) != 1:
        raise ValueError(f"{path}: expected exactly one statistics row for {sample_key!r}, got {len(found)}")
    return found[0]


def _featurecounts_summary(path, sample_key):
    header, rows = read_table(path)
    if header[0] != "Status" or sample_key not in header[1:]:
        raise ValueError(f"{path}: summary sample field {sample_key!r} absent")
    result = {}
    for row in rows:
        category = row["Status"]
        if not category or category in result:
            raise ValueError(f"{path}: duplicate/empty summary category")
        result[category] = number(row[sample_key], f"{path}: {category}", True)
    return result


def _metrics(input_row, sample, genes, special, stats_path, stats_key, stats_format):
    metric = dict(sample)
    metric.update(assay_id=input_row["assay_id"], source_project=input_row["source_project"],
                  quant_tool=input_row["format"].split("_")[0], qc_level="library",
                  gene_count_sum=total(genes.values()), genes_detected=sum(v > 0 for v in genes.values()),
                  quantified_gene_count=len(genes), gene_universe=input_row.get("gene_universe", "unknown"),
                  count_unit=sample.get("count_unit_override") or input_row.get("count_unit", "unknown"),
                  assignment_mode=input_row.get("assignment_mode", "unknown"),
                  assigned_count=None, counting_processed_count=None, assignment_rate_pct=None,
                  assignment_numerator=None, assignment_denominator=None, assignment_scope="counter_input",
                  assignment_counter_rate_pct=None, assignment_counter_status="unavailable",
                  assignment_counter_reason="No sufficient assignment statistics",
                  assignment_semantics_source="input declaration",
                  mapping_numerator=None, mapping_denominator=None, mapping_rate_pct=None,
                  mapping_scope=None, mapping_unit=None, mapping_status="unavailable",
                  mapping_reason="No mapping statistics provided", zero_assigned=None,
                  assignment_status="unavailable", assignment_reason="No sufficient assignment statistics",
                  stats_source=str(stats_path or ""), stats_sample_id=stats_key,
                  composition_valid=False)
    metric["genes_detected_pct"] = percentage(metric["genes_detected"], len(genes)) if metric["gene_universe"] == "full" else None
    metric["supplied_gene_detection_pct"] = percentage(metric["genes_detected"], len(genes))
    categories = dict(special)
    legacy = {}
    if stats_path and stats_format == "htseq_legacy":
        legacy = _legacy_row(stats_path, stats_key)
        for category in HTSEQ_COUNTERS:
            name = next((key for key in (category, "htseq_" + category, "__" + category) if key in legacy), None)
            if name and legacy[name] not in ("", "NA"):
                categories[category] = number(legacy[name], f"{stats_path}: {name}", True)
    elif stats_path and stats_format == "htseq_inline":
        _, categories = _htseq(stats_path, integer=input_row.get("assignment_mode") == "exclusive")
    elif stats_path and stats_format == "featurecounts_summary":
        categories = _featurecounts_summary(stats_path, stats_key)
    elif stats_path and stats_format == "salmon_meta":
        with Path(stats_path).open(encoding="utf-8") as handle:
            data = json.load(handle)
        if data.get("num_processed") is not None and data.get("num_mapped") is not None:
            mapped = number(data["num_mapped"], "Salmon num_mapped", True)
            processed = number(data["num_processed"], "Salmon num_processed", True)
            if mapped > processed:
                raise ValueError(f"{stats_path}: mapped fragments exceed processed fragments")
            metric.update(mapping_numerator=mapped, mapping_denominator=processed,
                          mapping_rate_pct=percentage(mapped, processed), mapping_scope="salmon_input",
                          mapping_unit="fragment", mapping_status="ok" if processed else "unavailable",
                          mapping_reason="" if processed else "Zero processed denominator")
    mode = input_row.get("assignment_mode", "unknown")
    complete_counts = metric["gene_universe"] == "full" or input_row["format"] in {"htseq_raw", "featurecounts_raw"}
    if metric["quant_tool"] == "htseq" and mode == "exclusive":
        reported_assigned = (number(legacy["assigned"], "legacy assigned")
                             if legacy.get("assigned") not in (None, "", "NA") else None)
        if complete_counts:
            metric["assigned_count"] = metric["gene_count_sum"]
            if reported_assigned is not None and reported_assigned != metric["gene_count_sum"]:
                raise ValueError(f"{stats_path}: assigned statistic disagrees with full count matrix for {stats_key}")
        else:
            metric["assigned_count"] = reported_assigned
        if all(category in categories for category in HTSEQ_COUNTERS) and metric["assigned_count"] is not None:
            metric["counting_processed_count"] = total([metric["assigned_count"], *categories.values()])
            if legacy.get("htseq_total") not in (None, "", "NA") and number(legacy["htseq_total"], "legacy HTSeq total") != metric["counting_processed_count"]:
                raise ValueError(f"{stats_path}: reported HTSeq total disagrees with counters for {stats_key}")
            metric["composition_valid"] = True
        else:
            reasons = []
            if metric["assigned_count"] is None:
                reasons.append("Gene universe is incomplete or unknown and no assigned-event statistic is supplied")
            if not all(category in categories for category in HTSEQ_COUNTERS):
                reasons.append("Incomplete HTSeq special counters; missing values are not zero")
            metric["assignment_reason"] = "; ".join(reasons)
    elif metric["quant_tool"] == "featurecounts" and "Assigned" in categories:
        metric["assigned_count"] = categories["Assigned"]
        if mode == "exclusive" and complete_counts and categories["Assigned"] != metric["gene_count_sum"]:
            raise ValueError(f"{stats_path}: featureCounts Assigned disagrees with exclusive count sum")
        if mode == "exclusive" and any(k.startswith("Unassigned_") for k in categories):
            metric["counting_processed_count"] = total(categories.values())
            metric["composition_valid"] = True
        else:
            metric["assignment_reason"] = "Summary counting semantics are not verified exclusive categories"
    elif mode != "exclusive":
        metric["assignment_reason"] = "Nonexclusive or unknown assignment policy: gene sum is not an assigned-event count"
    numerator, denominator = metric["assigned_count"], metric["counting_processed_count"]
    metric.update(assignment_numerator=numerator, assignment_denominator=denominator,
                  assignment_rate_pct=percentage(numerator, denominator),
                  zero_assigned=numerator == 0 if numerator is not None else None)
    if denominator is not None:
        if numerator > denominator:
            raise ValueError(f"{sample['sample_id']}: assigned count exceeds its denominator")
        metric["assignment_status"] = "ok" if denominator else "unavailable"
        metric["assignment_reason"] = "" if denominator else "Zero processed denominator"
    metric["assignment_counter_rate_pct"] = metric["assignment_rate_pct"]
    metric["assignment_counter_status"] = metric["assignment_status"]
    metric["assignment_counter_reason"] = ("Exclusive counting semantics declared; ratio uses each library's native counter events"
        if metric["assignment_counter_rate_pct"] is not None else metric["assignment_reason"])
    if metric["count_unit"] == "unknown":
        prior_reason = metric["assignment_reason"]
        metric.update(assignment_rate_pct=None, assignment_status="unavailable",
                      assignment_reason=(prior_reason + "; " if prior_reason else "") + "Unknown physical counting unit; native-counter fraction is reported separately when available")
        if metric["assignment_counter_rate_pct"] is not None:
            metric["assignment_counter_reason"] += "; read versus fragment unit is unverified; do not pool denominators across libraries"
    unassigned = [{"sample_id": sample["sample_id"], "assay_id": input_row["assay_id"],
                   "tool": metric["quant_tool"], "category": key, "count": value,
                   "unit": metric["count_unit"], "composition_valid": metric["composition_valid"],
                   "qc_level": "library"} for key, value in categories.items()]
    return metric, unassigned


def _technical_groups(sample_rows, all_ids):
    groups = defaultdict(list)
    for sample in sample_rows:
        groups[sample.get("technical_merge_id") or sample["sample_id"]].append(sample)
    excluded_keys = {"sample_id", "input_id", "source_sample_id", "stats_sample_id", "stats_source_path",
                     "technical_partition_id", "technical_disjoint", "source_project", "counts_source",
                     "run_metadata_path", "display_name", "count_unit_override"}
    for group_id, members in groups.items():
        if not members[0].get("technical_merge_id"):
            continue
        if group_id in all_ids:
            raise ValueError(f"Technical aggregate ID {group_id!r} collides with a raw sample ID")
        partitions = [s.get("technical_partition_id") for s in members]
        if any(not p for p in partitions) or len(set(partitions)) != len(partitions):
            raise ValueError(f"{group_id}: unique technical_partition_id values are required")
        if any(str(s.get("technical_disjoint", "")).lower() != "true" for s in members):
            raise ValueError(f"{group_id}: technical_disjoint=true is required")
        keys = set().union(*(s.keys() for s in members)) - excluded_keys
        for key in keys:
            if len({str(s.get(key, "")) for s in members}) > 1:
                raise ValueError(f"{group_id}: technical members disagree on {key}")
    return groups


def import_inputs(inputs, samples, base_dir, out_dir):
    """Import all raw samples and write assay matrices while preserving raw QC.

    An assay item has counts_path, samples, salmon (original files for tximport),
    errors, and gene_count. Only assays without errors have a counts_path.
    """
    base_dir, out_dir = Path(base_dir), Path(out_dir)
    result = {"assays": {}, "metrics": [], "unassigned": [], "provenance": [], "warnings": [], "errors": []}
    input_map = {row["input_id"]: row for row in inputs}
    if len(input_map) != len(inputs):
        raise ValueError("Duplicate input_id")
    ids = [s["sample_id"] for s in samples]
    if len(set(ids)) != len(ids):
        raise ValueError("Duplicate sample_id")
    if "gene_id" in ids:
        raise ValueError("sample_id gene_id conflicts with the count matrix identifier header")
    by_input = defaultdict(list)
    for row in samples:
        if row["input_id"] not in input_map:
            raise ValueError(f"Unknown input_id {row['input_id']!r}")
        by_input[row["input_id"]].append(dict(row))
    seen_paths, seen_stats, assay_counts, assay_inputs = set(), set(), defaultdict(dict), defaultdict(list)
    hash_cache = {}
    def file_hash(path):
        key = str(Path(path).resolve())
        if key not in hash_cache:
            hash_cache[key] = digest(path)
        return hash_cache[key]
    for input_row in inputs:
        input_id, assay = input_row["input_id"], input_row["assay_id"]
        input_row = dict(input_row)
        fmt = input_row["format"]
        if fmt not in FORMATS:
            raise ValueError(f"{input_id}: unsupported count format {fmt}")
        allowed_stats = {"htseq": {"none", "htseq_inline", "htseq_legacy"},
                         "featurecounts": {"none", "featurecounts_summary"},
                         "salmon": {"none", "salmon_meta"}, "generic": {"none"}}
        if input_row.get("stats_format", "none") not in allowed_stats[fmt.split("_")[0]]:
            raise ValueError(f"{input_id}: stats_format is incompatible with {fmt}")
        path = resolve(base_dir, input_row["counts_path"])
        canonical = str(path.resolve())
        if canonical in seen_paths:
            raise ValueError(f"Count file registered more than once: {path}")
        seen_paths.add(canonical)
        members = by_input[input_id]
        if not members:
            raise ValueError(f"{input_id}: no sample mappings")
        integer = input_row.get("assignment_mode") == "exclusive" and fmt != "salmon_quant"
        special, salmon_rows = {}, []
        if fmt == "htseq_raw":
            if len(members) != 1:
                raise ValueError(f"{input_id}: HTSeq raw file requires exactly one sample mapping")
            genes, special = _htseq(path, integer)
            matrix = {members[0]["source_sample_id"]: genes}
        elif fmt == "salmon_quant":
            if len(members) != 1:
                raise ValueError(f"{input_id}: Salmon quant.sf requires one sample mapping")
            map_path = resolve(base_dir, input_row["gene_map_path"]) if input_row.get("gene_map_path") else None
            genes = _salmon(path, map_path)
            matrix = {members[0]["source_sample_id"]: genes}
            salmon_rows = [{"sample_id": members[0]["sample_id"], "quant_path": str(path.resolve()), "gene_map_path": str(map_path.resolve())}]
        else:
            matrix = _matrix(path, fmt, integer)
        expected = [s["source_sample_id"] for s in members]
        if len(set(expected)) != len(expected) or set(expected) != set(matrix):
            raise ValueError(f"{input_id}: count/sample mapping is not 1:1; unmapped columns={sorted(set(matrix)-set(expected))}, absent columns={sorted(set(expected)-set(matrix))}")
        item = result["assays"].setdefault(assay, {"counts_path": "", "samples": [], "salmon": [], "errors": []})
        item["salmon"].extend(salmon_rows)
        provenance = dict(input_row, counts_path=canonical, counts_sha256=file_hash(path))
        provenance["effective_provenance_level"] = "declared"
        if input_row.get("run_metadata_path"):
            meta_path = resolve(base_dir, input_row["run_metadata_path"])
            provenance["run_metadata_sha256"] = file_hash(meta_path)
            provenance["run_metadata_path"] = str(meta_path.resolve())
            if meta_path.suffix.lower() == ".json":
                with meta_path.open(encoding="utf-8-sig") as handle:
                    metadata = json.load(handle)
                if not isinstance(metadata, dict):
                    raise ValueError(f"{meta_path}: run metadata must be a JSON object")
                options = metadata.get("quantification", metadata.get("options", metadata))
                if isinstance(options, dict):
                    for key in ("count_feature", "count_attribute", "nonunique", "count_multimappers", "overlap_mode", "fractional", "tool_version"):
                        if key in options:
                            input_row["effective_" + key] = json.dumps(options[key], sort_keys=True)
                            provenance["effective_" + key] = options[key]
                    for key in ("assignment_mode", "count_unit"):
                        if options.get(key) and input_row.get(key) not in ("", "unknown", options[key]):
                            raise ValueError(f"{input_id}: declared {key} conflicts with run metadata")
                    if options.get("nonunique") in ("all", "fraction", "random") and input_row.get("assignment_mode") == "exclusive":
                        raise ValueError(f"{input_id}: exclusive assignment conflicts with HTSeq nonunique option")
                provenance["run_metadata_json"] = json.dumps(metadata, ensure_ascii=False, sort_keys=True)
        if input_row.get("gene_map_path"):
            provenance["gene_map_sha256"] = file_hash(resolve(base_dir, input_row["gene_map_path"]))
        assay_inputs[assay].append(input_row)
        result["provenance"].append(provenance)
        if input_row.get("provenance_level") == "verified":
            result["warnings"].append(f"{input_id}: user verified declaration retained; reference content not independently checked")
        for sample in members:
            sample.update(assay_id=assay, source_project=input_row["source_project"],
                          counts_source=canonical, count_unit=sample.get("count_unit_override") or input_row.get("count_unit", "unknown"))
            genes = matrix[sample["source_sample_id"]]
            stats_value = sample.get("stats_source_path") or input_row.get("stats_path")
            stats_path = resolve(base_dir, stats_value) if stats_value else None
            stats_key = sample.get("stats_sample_id") or sample["source_sample_id"]
            stats_format = input_row.get("stats_format", "none")
            if stats_path:
                stats_join = (str(stats_path.resolve()), stats_key)
                if stats_join in seen_stats:
                    raise ValueError(f"Statistics row is reused by multiple samples: {stats_join}")
                seen_stats.add(stats_join)
                result["provenance"].append({"sample_id": sample["sample_id"], "stats_path": str(stats_path.resolve()), "stats_sample_id": stats_key, "stats_sha256": file_hash(stats_path)})
            metric, categories = _metrics(input_row, sample, genes, special, stats_path, stats_key, stats_format)
            result["metrics"].append(metric)
            result["unassigned"].extend(categories)
            assay_counts[assay][sample["sample_id"]] = genes
            item["samples"].append(sample)
    metric_lookup = {m["sample_id"]: m for m in result["metrics"]}
    for assay, item in result["assays"].items():
        entries, per_sample = assay_inputs[assay], assay_counts[assay]
        for key in ("reference_id", "annotation_id", "assignment_mode", "gene_universe"):
            if len({str(entry.get(key, "")) for entry in entries}) > 1:
                item["errors"].append(f"{assay}: incompatible {key}")
        for key in ("reference_sha256", "annotation_sha256"):
            known = {entry[key] for entry in entries if entry.get(key)}
            if len(known) > 1:
                item["errors"].append(f"{assay}: conflicting {key}")
        for key in ("count_feature", "count_attribute", "nonunique", "count_multimappers", "overlap_mode", "fractional"):
            known = {entry["effective_" + key] for entry in entries if "effective_" + key in entry}
            if len(known) > 1:
                item["errors"].append(f"{assay}: conflicting quantification option {key}")
        tools = {entry["format"].split("_")[0] for entry in entries}
        if len(tools) != 1:
            item["errors"].append(f"{assay}: different quantification tools cannot share an assay")
        if any(isinstance(value, int) and value > 2**53 for genes in per_sample.values() for value in genes.values()):
            item["errors"].append(f"{assay}: count exceeds exact double-precision integer range for R analysis; exact QC retained")
        gene_order = list(next(iter(per_sample.values())))
        if any(set(genes) != set(gene_order) for genes in per_sample.values()):
            item["errors"].append(f"{assay}: gene universes differ; intersection/zero-fill is forbidden")
        groups = _technical_groups(item["samples"], set(ids))
        if item["errors"]:
            result["errors"].extend(item["errors"])
            continue
        aggregate_counts, aggregate_samples = {}, []
        for group_id, members in groups.items():
            members_ids = [s["sample_id"] for s in members]
            combined = {gene: total(per_sample[sid][gene] for sid in members_ids) for gene in gene_order}
            merged = dict(members[0], sample_id=group_id, constituent_sample_ids="|".join(members_ids))
            if members[0].get("technical_merge_id"):
                if "salmon" in tools:
                    item["errors"].append(f"{assay}: Salmon technical aggregation requires an explicit tximport-compatible workflow")
                    continue
                combined_units = {metric_lookup[sid]["count_unit"] for sid in members_ids}
                if len(combined_units) != 1:
                    raise ValueError(f"{group_id}: incompatible technical count units")
                raw_metrics = [metric_lookup[sid] for sid in members_ids]
                metric = dict(raw_metrics[0], **{k: v for k, v in merged.items() if k not in raw_metrics[0]})
                metric.update(sample_id=group_id, qc_level="technical_aggregate", constituent_sample_ids="|".join(members_ids),
                              gene_count_sum=total(combined.values()), genes_detected=sum(v > 0 for v in combined.values()))
                metric["genes_detected_pct"] = percentage(metric["genes_detected"], len(gene_order)) if metric["gene_universe"] == "full" else None
                metric["supplied_gene_detection_pct"] = percentage(metric["genes_detected"], len(gene_order))
                for key in ("assigned_count", "counting_processed_count", "assignment_numerator", "assignment_denominator", "mapping_numerator", "mapping_denominator"):
                    metric[key] = total(m[key] for m in raw_metrics) if all(m.get(key) is not None for m in raw_metrics) else None
                metric["assignment_rate_pct"] = (percentage(metric["assignment_numerator"], metric["assignment_denominator"])
                                                  if all(m["assignment_status"] == "ok" for m in raw_metrics) else None)
                metric["mapping_rate_pct"] = (percentage(metric["mapping_numerator"], metric["mapping_denominator"])
                                               if all(m["mapping_status"] == "ok" for m in raw_metrics) else None)
                metric["mapping_status"] = "ok" if metric["mapping_rate_pct"] is not None else "unavailable"
                metric["mapping_reason"] = "" if metric["mapping_status"] == "ok" else "One or more technical components lack mapping statistics"
                metric["zero_assigned"] = metric["assigned_count"] == 0 if metric["assigned_count"] is not None else None
                metric["assignment_status"] = "ok" if metric["assignment_rate_pct"] is not None else "unavailable"
                metric["assignment_reason"] = "" if metric["assignment_status"] == "ok" else "One or more technical components lack a compatible denominator"
                counter_compatible = (metric["count_unit"] != "unknown" and
                                      all(m["assignment_counter_status"] == "ok" for m in raw_metrics))
                metric["assignment_counter_rate_pct"] = (percentage(metric["assignment_numerator"], metric["assignment_denominator"])
                                                         if counter_compatible else None)
                metric["assignment_counter_status"] = "ok" if metric["assignment_counter_rate_pct"] is not None else "unavailable"
                metric["assignment_counter_reason"] = ("Compatible known units; exclusive counter counts summed before division"
                    if metric["assignment_counter_status"] == "ok" else "Technical aggregate needs known compatible counting units and complete component counters")
                metric["composition_valid"] = all(m["composition_valid"] for m in raw_metrics)
                metric["source_project"] = "|".join(dict.fromkeys(m["source_project"] for m in members))
                metric["input_id"] = "|".join(dict.fromkeys(m["input_id"] for m in members))
                metric["stats_source"] = "|".join(dict.fromkeys(m["stats_source"] for m in raw_metrics))
                result["metrics"].append(metric)
                raw_categories = [r for r in result["unassigned"] if r["sample_id"] in members_ids and r["qc_level"] == "library"]
                for category in dict.fromkeys(row["category"] for row in raw_categories):
                    parts = [row for row in raw_categories if row["category"] == category]
                    result["unassigned"].append(dict(sample_id=group_id, assay_id=assay, tool=metric["quant_tool"],
                        category=category, count=total(r["count"] for r in parts) if len(parts) == len(members_ids) else None,
                        unit=metric["count_unit"], composition_valid=metric["composition_valid"], qc_level="technical_aggregate"))
                merged["source_project"] = "|".join(dict.fromkeys(m["source_project"] for m in members))
            aggregate_counts[group_id] = combined
            aggregate_samples.append(merged)
        if any(isinstance(value, int) and value > 2**53 for values in aggregate_counts.values() for value in values.values()):
            item["errors"].append(f"{assay}: technical aggregate exceeds exact double-precision integer range; QC retained")
        item["samples"] = aggregate_samples
        item["gene_count"] = len(gene_order)
        if item["errors"]:
            result["errors"].extend(item["errors"])
            continue
        matrix_path = out_dir / "imported" / assay / "counts.tsv"
        item["counts_path"] = write_tsv(matrix_path,
            (dict(gene_id=gene, **{sid: values[gene] for sid, values in aggregate_counts.items()}) for gene in gene_order),
            ["gene_id", *aggregate_counts])
        item["samples_path"] = write_tsv(matrix_path.parent / "samples.tsv", aggregate_samples)
    return result
