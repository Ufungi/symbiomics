"""Write auditable QC tables; visualization/report generation is separate."""
from pathlib import Path

try:
    from .importers import write_tsv
except ImportError:
    from importers import write_tsv

DEFINITIONS = [
    dict(metric="gene_count_sum", numerator="sum of unfiltered supplied gene counts", denominator="", unit="count_unit", scope="supplied count matrix"),
    dict(metric="assigned_count", numerator="assigned events from valid tool semantics", denominator="", unit="count_unit", scope="counter input"),
    dict(metric="assignment_rate_pct", numerator="assignment_numerator", denominator="assignment_denominator", unit="percent; event unit in count_unit", scope="assignment_scope; never assumed sequenced input"),
    dict(metric="assignment_counter_rate_pct", numerator="assignment_numerator", denominator="assignment_denominator", unit="percent of native counter events", scope="Within-library exclusive-counter declaration; physical read/fragment unit may be unknown; never pool unknown units"),
    dict(metric="genes_detected", numerator="number of supplied genes with count > 0", denominator="", unit="gene", scope="before filtering"),
    dict(metric="genes_detected_pct", numerator="genes_detected", denominator="quantified_gene_count", unit="percent", scope="full quantified universe only; otherwise NA"),
    dict(metric="supplied_gene_detection_pct", numerator="genes_detected", denominator="quantified_gene_count", unit="percent", scope="supplied matrix rows; may be prefiltered"),
    dict(metric="mapping_rate_pct", numerator="mapping_numerator", denominator="mapping_denominator", unit="percent; event unit in mapping_unit", scope="mapping_scope"),
    dict(metric="zero_assigned", numerator="assigned_count == 0", denominator="", unit="boolean", scope="actual unrounded count; NA if unknown"),
]


def emit_qc(result, out_dir, plots=None):
    """Write raw/aggregated metrics and return absolute paths by artifact name.

    out_dir is the run directory; outputs go in 00_quant_qc. Every rate retains
    numeric numerator/denominator, unit, scope, source and status/reason columns.
    plots is accepted for API compatibility; report code owns visualization.
    """
    root = Path(out_dir) / "00_quant_qc"
    paths = {
        "sample_metrics": write_tsv(root / "sample_metrics.tsv", result.get("metrics", [])),
        "metric_definitions": write_tsv(root / "metric_definitions.tsv", DEFINITIONS),
        "unassigned": write_tsv(root / "unassigned_long.tsv", result.get("unassigned", []),
                                ["sample_id", "assay_id", "tool", "category", "count", "unit", "composition_valid", "qc_level"]),
        "provenance": write_tsv(root / "provenance.tsv", result.get("provenance", [])),
    }
    return paths
