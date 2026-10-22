"""Validated table contracts. Data cells are never evaluated as code."""
from __future__ import annotations
import csv
import math
import re
from pathlib import Path

ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
COLUMN = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
INPUT_COLUMNS = """input_id source_project assay_id format counts_path stats_path stats_format gene_universe run_metadata_path reference_id annotation_id reference_sha256 annotation_sha256 count_unit assignment_mode provenance_level gene_map_path""".split()
SAMPLE_COLUMNS = """sample_id input_id source_sample_id stats_sample_id count_unit_override stats_source_path biological_unit_id library_id technical_merge_id technical_partition_id technical_disjoint study_id sample_sets include_analysis exclusion_reason subject_id""".split()
ANALYSIS_COLUMNS = """analysis_id sample_set assay_id analysis_mode design_type factor_cols numeric_cols nuisance_factor_cols nuisance_numeric_cols interaction_terms subject_col reference_levels numeric_units normalization fdr_cut lfc_cut effect_test enrichment go_method pathway_methods adjusted_view adjusted_batch_col seed""".split()
CONTRAST_COLUMNS = "analysis_id contrast_id basis term weight test_id".split()
ANNOTATION_COLUMNS = "annotation_id assay_id kind path id_map_path version source universe_path".split()
PLOT_COLUMNS = "plot_id scope analysis_id metric x color facet order scale".split()
COLUMNS = dict(inputs=INPUT_COLUMNS, samples=SAMPLE_COLUMNS, analyses=ANALYSIS_COLUMNS,
               contrasts=CONTRAST_COLUMNS, annotations=ANNOTATION_COLUMNS, plots=PLOT_COLUMNS)
ANALYSIS_DEFAULTS = dict(analysis_mode="inferential", normalization="TMM", fdr_cut="0.05",
                        lfc_cut="1", effect_test="qlf", enrichment="none",
                        adjusted_view="none", seed="123")
FORMATS = {"htseq_matrix", "featurecounts_matrix", "htseq_raw", "featurecounts_raw",
           "salmon_quant", "generic_raw_matrix"}

def safe_id(value, field="id"):
    if not ID.fullmatch(value or "") or ".." in value:
        raise ValueError(f"{field}: unsafe or empty identifier {value!r}")
    return value

def bool_value(value, field):
    if str(value).lower() not in {"true", "false"}:
        raise ValueError(f"{field}: expected true/false, got {value!r}")
    return str(value).lower() == "true"

def path_value(value, base):
    """Accept Windows drive paths when invoked in WSL; never invoke a shell."""
    s = str(value)
    if re.match(r"^[A-Za-z]:[\\/]", s) and __import__("os").name != "nt":
        s = "/mnt/" + s[0].lower() + "/" + s[3:].replace(chr(92), "/")
    return (Path(s) if Path(s).is_absolute() else Path(base) / s).resolve()

def read_table(path, name, optional=False):
    path = Path(path)
    if not path.exists():
        if optional:
            return []
        raise ValueError(f"missing required table: {path}")
    with path.open(encoding="utf-8-sig", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError(f"empty table without header: {path}")
        if any(not k for k in header) or len(set(header)) != len(header):
            raise ValueError(f"{path}: blank or duplicate header")
        unknown = set(header) - set(COLUMNS[name])
        if name != "samples" and unknown:
            raise ValueError(f"{path}: unknown columns: {sorted(unknown)}")
        if name == "samples" and any(not COLUMN.fullmatch(k) for k in header):
            raise ValueError(f"{path}: metadata column names must be safe identifiers")
        rows = []
        for line, cells in enumerate(reader, 2):
            if not cells or not any(cells):
                continue
            if len(cells) != len(header):
                raise ValueError(f"{path}:{line}: expected {len(header)} cells, got {len(cells)}")
            row = dict(zip(header, cells))
            if any("\x00" in v for v in row.values()):
                raise ValueError(f"{path}:{line}: NUL in input")
            rows.append(row)
    return rows

def required(row, fields, label):
    missing = [k for k in fields if not row.get(k, "")]
    if missing:
        raise ValueError(f"{label}: missing {', '.join(missing)}")

def validate_import_tables(inputs, samples):
    if not inputs or not samples:
        raise ValueError("inputs.tsv and samples.tsv must contain at least one data row")
    seen = set()
    for row in inputs:
        required(row, ["input_id","source_project","assay_id","format","counts_path",
                       "reference_id","annotation_id","count_unit","assignment_mode"], "input")
        for key in ("input_id", "assay_id"):
            safe_id(row[key], key)
        if row["input_id"] in seen:
            raise ValueError(f"duplicate input_id: {row['input_id']}")
        seen.add(row["input_id"])
        for hash_key in ("reference_sha256","annotation_sha256"):
            if row.get(hash_key) and not re.fullmatch(r"[0-9A-Fa-f]{64}",row[hash_key]):
                raise ValueError(f"invalid SHA256 in {hash_key}")
        if row["format"] not in FORMATS:
            raise ValueError(f"unsupported format: {row['format']}")
        for key, allowed in {
            "count_unit": {"read","fragment","estimated_fragment","unknown"},
            "assignment_mode": {"exclusive","multiassign","fractional","unknown"},
            "gene_universe": {"full","prefiltered","unknown"},
            "stats_format": {"htseq_inline","htseq_legacy","featurecounts_summary","salmon_meta","none"},
            "provenance_level": {"verified","declared"},
        }.items():
            if key not in row or row[key] == "":
                row[key] = {"gene_universe":"unknown","stats_format":"none","provenance_level":"declared"}.get(key,"")
            if row[key] not in allowed:
                raise ValueError(f"input {row['input_id']}: invalid {key}={row[key]!r}")
    ids, links = set(), set()
    input_assays={r["input_id"]:r["assay_id"] for r in inputs}
    merge_assays={}
    for row in samples:
        required(row, ["sample_id","input_id","source_sample_id","biological_unit_id",
                       "library_id","study_id","include_analysis"], "sample")
        safe_id(row["sample_id"], "sample_id")
        if row["sample_id"] in ids:
            raise ValueError(f"duplicate sample_id: {row['sample_id']}")
        ids.add(row["sample_id"])
        if row["input_id"] not in seen:
            raise ValueError(f"unknown input_id: {row['input_id']}")
        link = (row["input_id"], row["source_sample_id"])
        if link in links:
            raise ValueError(f"count column mapped twice: {link}")
        links.add(link)
        if not bool_value(row["include_analysis"], "include_analysis") and not row.get("exclusion_reason"):
            raise ValueError(f"sample {row['sample_id']}: excluded sample needs exclusion_reason")
        for group in filter(None, row.get("sample_sets","").split("|")):
            safe_id(group, "sample_set")
        if row.get("count_unit_override") and row["count_unit_override"] not in {"read","fragment","estimated_fragment","unknown"}:
            raise ValueError("invalid count_unit_override")
        if row.get("biological_unit_status","") not in {"","verified","unverified","unknown"}:
            raise ValueError("biological_unit_status must be verified/unverified/unknown")
        if row.get("technical_merge_id"):
            mid=row["technical_merge_id"];assay=input_assays[row["input_id"]]
            if mid in merge_assays and merge_assays[mid]!=assay:raise ValueError("technical_merge_id reused across assays")
            merge_assays[mid]=assay
            safe_id(row["technical_merge_id"], "technical_merge_id")
            required(row, ["technical_partition_id","technical_disjoint"], "technical aggregation")
            if not bool_value(row["technical_disjoint"], "technical_disjoint"):
                raise ValueError("technical aggregation requires disjoint data declaration")

    if set(merge_assays) & ids:
        raise ValueError("technical aggregate ID collides with raw sample ID")

def analysis_defaults(row):
    out = {**{k:"" for k in ANALYSIS_COLUMNS}, **row}
    for key, value in ANALYSIS_DEFAULTS.items():
        if not out.get(key):
            out[key] = value
    return out

def validate_analysis(row, contrasts, samples, allow_empty_contrasts=False):
    required(row, ["analysis_id","sample_set","assay_id","design_type"], "analysis")
    for key in ("analysis_id","sample_set","assay_id"):
        safe_id(row[key],key)
    if row["analysis_mode"] not in {"inferential","exploratory"}:
        raise ValueError("analysis_mode must be inferential or exploratory")
    if row["design_type"] not in {"group","additive","factorial","paired","trend"}:
        raise ValueError("unsupported design_type")
    if row["normalization"] != "TMM" or row["effect_test"] not in {"qlf","treat","deseq2"}:
        raise ValueError("normalization=TMM and effect_test=qlf|treat|deseq2 are supported")
    for key, lower, upper in [("fdr_cut",0,1),("lfc_cut",0,float("inf"))]:
        try: value = float(row[key])
        except ValueError: raise ValueError(f"invalid {key}")
        if not math.isfinite(value) or value < lower or value > upper or (key=="fdr_cut" and value==0):
            raise ValueError(f"invalid {key}")
    try:
        seed=int(row["seed"])
        if not 0 <= seed <= 2147483647: raise ValueError()
    except ValueError: raise ValueError("seed must be an integer in 0..2147483647")
    declared=[]
    for key in ("factor_cols","numeric_cols","nuisance_factor_cols","nuisance_numeric_cols"):
        vals = list(filter(None,row.get(key,"").split("|")))
        if len(vals)!=len(set(vals)): raise ValueError(f"duplicate variable in {key}")
        declared += vals
    if row.get("subject_col"): declared.append(row["subject_col"])
    if len(declared)!=len(set(declared)):
        raise ValueError("a variable was declared more than once")
    if any(not COLUMN.fullmatch(k) for k in declared):
        raise ValueError("invalid model variable name")
    for key in declared:
        if any(key not in s or s[key] == "" for s in samples):
            raise ValueError(f"model metadata missing: {key}")
    if row["analysis_mode"] == "inferential" and not contrasts and not allow_empty_contrasts:
        raise ValueError("inferential analysis needs declared contrasts")
    if row["analysis_mode"] == "exploratory" and contrasts:
        raise ValueError("exploratory analysis must not declare inferential contrasts")
    if row["adjusted_view"] not in {"none","combat_seq"}:
        raise ValueError("unsupported adjusted_view")
    if row["adjusted_view"]=="combat_seq" and not row["adjusted_batch_col"]:
        raise ValueError("combat_seq needs adjusted_batch_col")
    if row["enrichment"] != "none":
        kinds=set(row["enrichment"].split("|"))
        if not kinds <= {"go","pathway"}: raise ValueError("unsupported enrichment")
        if "go" in kinds and row["go_method"] != "fisher_ora":
            raise ValueError("GO enrichment needs go_method=fisher_ora")
        if "pathway" in kinds and (not row["pathway_methods"] or
                not set(row["pathway_methods"].split("|")) <= {"camera","fgsea"}):
            raise ValueError("pathway_methods must explicitly name camera/fgsea")
    groups={}
    for c in contrasts:
        required(c, ["contrast_id","basis","term","weight"], "contrast")
        safe_id(c["contrast_id"],"contrast_id")
        if c.get("test_id"): safe_id(c["test_id"],"test_id")
        if c["basis"] not in {"cell","coefficient"}: raise ValueError("invalid contrast basis")
        try: w=float(c["weight"])
        except ValueError: raise ValueError("contrast weight is not numeric")
        if not math.isfinite(w): raise ValueError("contrast weight must be finite")
        groups.setdefault(c["contrast_id"],[]).append(c)
    for cid, rows in groups.items():
        if len({(c["basis"],c.get("test_id","")) for c in rows}) != 1:
            raise ValueError(f"{cid}: inconsistent basis/test_id")
        if rows[0]["basis"]=="cell" and abs(sum(float(c["weight"]) for c in rows))>1e-10:
            raise ValueError(f"{cid}: cell weights must sum to zero")
        if row["effect_test"]=="deseq2" and any(c.get("test_id") for c in rows):
            raise ValueError("effect_test=deseq2 does not support omnibus test_id rows")
        if row["effect_test"]=="treat" and any(c.get("test_id") for c in rows):
            raise ValueError("omnibus+treat is unsupported")
