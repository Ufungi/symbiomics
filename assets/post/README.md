# Post-analysis project

Replace all example rows with your own inputs and metadata; no example counts are supplied.
1. inputs.tsv: raw gene count files and optional original quantification summaries.
2. samples.tsv: exact count columns, metadata and biological units. Never infer replicates from file count.
3. analyses.tsv: independent, paired, factorial or numeric models. Defaults are only examples.
4. contrasts.tsv: explicit model comparisons; positive weights define the positive effect.
5. annotations.tsv / plots.tsv: optional annotation sources and figure choices.

Run `symbiomics post validate --project-dir .` to inspect model matrices and contrasts.
Run `symbiomics post run --project-dir .` for QC, models and reports.
Run `symbiomics post run --project-dir . --qc-only` when only quantification QC is wanted.
For PCA without inferential testing set analysis_mode=exploratory and remove its contrast rows.

Paths are relative to each table. Raw QC retains excluded and all-zero samples.
QC assignment requires original summaries and compatible denominator units; missing metrics remain NA.
Outcomes appear in output/runs/<run_id>/report.html; previous runs are preserved.
Full documentation: docs/post_analysis.md in the Symbiomics repository.
