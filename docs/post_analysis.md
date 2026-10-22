# Table-driven RNA-seq post-analysis

`symbiomics post` consumes existing gene counts and experimental metadata. It
writes quantification QC before statistical model validation, followed by
exploratory plots and, when a design is confirmed, edgeR analysis. It is a
standalone command family: creating tables and running quantification QC do
not start Java, Nextflow, alignment or annotation jobs.

The tables describe the experiment. Sample names are never used to guess
control groups, tissue, treatment, pairing or independent biological replicates.
The same commands work with separate execution projects when their count
semantics and reference/annotation are compatible.

## Start with quantification QC

Run from the repository or use an installed `symbiomics` launcher:

```bash
scripts/symbiomics post init -p my_study
# Edit post/my_study/inputs.tsv and samples.tsv.
scripts/symbiomics post validate -p my_study --qc-only
scripts/symbiomics post run -p my_study --qc-only
```

`-p NAME` selects `post/NAME/` inside this repository. `--project-dir PATH`
selects an explicit project directory; the two selectors are alternatives.
`post init` refuses to overwrite an existing project.

```bash
scripts/symbiomics post init --project-dir /data/experiments/fungal_qc
scripts/symbiomics post run --project-dir /data/experiments/fungal_qc --qc-only
```

For a Windows data file used from WSL, either use its `/mnt/c/...` path or a
Windows drive path in the input table. For example, the counting matrix
`C:\AI\Research\Piden_transcriptome\htseq-stat\piden2trima_all_htseq.csv`
can be entered as
`/mnt/c/AI/Research/Piden_transcriptome/htseq-stat/piden2trima_all_htseq.csv`.
The matrix format and counting units must still be declared from its actual
provenance. A `.csv` extension alone does not establish counting semantics.

Open the `report.html` path printed by the command. Its first section includes
all imported observations, including zero libraries and samples excluded from
later analysis. No infection, colonization or general sample-quality threshold
is inferred from low assignment alone.

## Fill the six standard tables

TSV files are UTF-8, with one header row and literal tab delimiters. Keep IDs
unique. Blank cells mean unavailable; do not replace unavailable counts with
zero. Paths are relative to the directory containing the project tables unless
absolute. Preserve gene and sample ID case exactly.

| Table | One row describes | Important fields |
| --- | --- | --- |
| `inputs.tsv` | One count input | `input_id`, `source_project`, `assay_id`, `format`, `counts_path`, `stats_path`, `stats_format`, `reference_id`, `annotation_id`, `count_unit`, `assignment_mode`, `gene_universe` |
| `samples.tsv` | One library or technical partition | `sample_id`, `input_id`, `source_sample_id`, `biological_unit_id`, `library_id`, `study_id`, `sample_sets`, `include_analysis`; arbitrary metadata columns such as `treatment`, `tissue`, `batch`, `time_h` |
| `analyses.tsv` | One model or exploratory selection | `analysis_id`, `sample_set`, `assay_id`, `analysis_mode`, `design_type`, factor/numeric/nuisance columns, reference levels, normalization and inference settings |
| `contrasts.tsv` | One weighted term of a comparison | `analysis_id`, `contrast_id`, `basis`, `term`, `weight`; optional `test_id` for a joint test |
| `annotations.tsv` | One optional gene annotation mapping | `annotation_id`, `assay_id`, `kind`, `path`, optional ID map and background universe |
| `plots.tsv` | One optional chart configuration | `plot_id`, `scope`, `analysis_id`, `metric`, `x`, `color`, `facet`, `order`, `scale` |

`inputs.tsv` and `samples.tsv` are required. Quantification QC can run with no
analysis or contrast rows. Header-only optional tables are valid. The scaffold
contains the complete column sets and a short example; replace example paths
and metadata with actual data before running.

Every count column must map exactly once through
`(input_id, source_sample_id)`. Map differently named summary columns explicitly
with `stats_sample_id`. `include_analysis=false` requires an
`exclusion_reason`; that sample remains in QC. If a matrix has already been
filtered, declare `gene_universe=prefiltered`. Use `unknown` when its completeness
is unconfirmed. Only `full` supports the full quantified-gene detection rate.

`source_project` records the execution source, `study_id` identifies a research
study, and `sample_sets` selects analyses. These fields need not be identical.
Separate sample-set IDs with `|`. A biological unit is an independent organism,
culture or experimental unit; different tissues from one organism are not
independent units simply because they have different sample IDs.

## Quantification inputs and rate definitions

| Format | Input | Quantification QC |
| --- | --- | --- |
| `htseq_matrix` | Gene rows and sample columns; provide counters in a separate statistics file | Counts and gene detection; assignment rate when exclusive semantics, known unit and complete special counters support a denominator |
| `htseq_raw` | Original two-column HTSeq output | Same rules, preserving original special counters |
| `featurecounts_matrix` | Merged gene count matrix | Gene counts/detection; use the original summary for assignment statistics |
| `featurecounts_raw` | Original featureCounts output | Gene counts plus supplied `.summary` categories |
| `generic_raw_matrix` | Nonnegative raw gene count matrix | Gene counts/detection; no assignment rate inferred from the matrix alone |
| `salmon_quant` | Original `quant.sf` plus transcript-to-gene mapping | Mapping statistics from `meta_info.json`; gene analysis uses the dedicated tximport route |

`stats_format` identifies `htseq_inline`, `htseq_legacy`,
`featurecounts_summary`, `salmon_meta` or `none`. Separate sample-specific
statistics can be supplied with `stats_source_path` in `samples.tsv`. Matrix inputs reject `__*` counter rows;
original two-column HTSeq files use `htseq_raw`, or can serve as a separate
`stats_path` with `stats_format=htseq_inline`.

For prefiltered or unknown-completeness matrices, an assignment numerator
must come from an independently reported assigned statistic (for example a
legacy `assigned` field or featureCounts summary). The supplied gene sum does
not reconstruct that numerator. Full-universe gene detection uses the explicit
`gene_universe=full` completeness declaration; other inputs report detection
within supplied rows separately.

* **Gene assignment rate:** assigned units divided by units processed at the
  counting stage, multiplied by 100. This is not automatically a fraction of
  the original sequenced reads. Unknown units, missing counters or
  incompatible assignment semantics produce `NA` and a reason.
* **Native counter assignment (`assignment_counter_rate_pct`):** assigned
  counters divided by processed counters within one library, multiplied by
  100, when declared exclusive/compatible counting semantics support that
  fraction. This can be available even when the original physical unit
  (read versus fragment) is unverified; `assignment_rate_pct` stays `NA` in
  that case. Status/reason fields retain the assumptions. This is a
  counting-stage proportion; a sequencing read mapping rate requires its own
  mapping-log denominator. Unknown-unit technical aggregates keep this metric
  unavailable, and reports do not pool unknown-unit count totals.
* **Detected-gene rate:** genes with count greater than zero divided by the
  full quantified gene universe. Detection within prefiltered or uncertain
  supplied rows is shown separately as `supplied_gene_detection_pct`.
* **Alignment/mapping rate:** mapped divided by processed units from mapping
  logs. It is distinct from gene assignment.
* **Gene count sum:** the sum of supplied counts. Multiassignment, fractional
  assignment and estimated counts can make this differ from the number of
  uniquely assigned reads or fragments.

Use `read`, `fragment`, `estimated_fragment` or `unknown` for `count_unit`, and
`exclusive`, `multiassign`, `fractional` or `unknown` for `assignment_mode`.
Do not choose `exclusive` merely to obtain a nonmissing rate. The report keeps
assay, counting unit, processing scope, tool, gene universe and raw/technical
aggregate QC separate. It preserves small positive values in scientific
notation and places true zero in a separate lane on log plots. Missing values
are never plotted as zero.

Counter category charts show original counts. Categories are not silently
completed with zeros or treated as an exhaustive, disjoint 100% composition.
Group summaries are descriptive: valid/missing counts, actual zeros,
quartiles, extrema and declared biological-unit counts. No automatic
ANOVA/Tukey tests are added.

## Metadata-based QC plots

For example, an assignment-rate chart grouped by tissue and colored by a
metadata field `condition`, split by source project:

```tsv
plot_id	scope	analysis_id	metric	x	color	facet	order	scale
assignment_by_tissue	quant_qc		assignment_rate_pct	tissue	condition	source_project	root|stem|needle	log10
```

Use actual tab characters in the file. Supported QC metrics include
`assignment_rate_pct`, `assignment_counter_rate_pct`, `genes_detected_pct`, `assigned_count`,
`gene_count_sum`, `genes_detected`, `quantified_gene_count`,
`supplied_gene_detection_pct` and `mapping_rate_pct`. Scales are `linear` or
`log10`. Each SVG has a TSV sidecar containing the sample values used to draw
it. Missing grouping metadata is labeled explicitly. The HTML includes inline
summary tables and local artifact links; it loads no network resources.

For analysis-coordinate charts, use `scope=analysis`, the selected
`analysis_id`, and `metric=pca` or `metric=mds`:

```tsv
plot_id	scope	analysis_id	metric	x	color	facet	order	scale
pca_by_tissue	analysis	main	pca		tissue	source_project	root|stem|needle	linear
mds_by_condition	analysis	main	mds		condition		control|treated	linear
```

These projections use the engine's saved coordinates and join the recorded
analysis sample metadata by `sample_id`; they never recalculate PCA/MDS from
QC values. Axes remain PC1/PC2 or MDS1/MDS2, so `x` is ignored and should be
left blank. `scale` is `linear` (or blank). `color` selects a metadata column,
`facet` splits panels, and `order` controls color/legend levels. Facets retain
the same coordinate limits and colors. Missing labels are explicit. Unknown
metadata columns, mismatched IDs and unsupported metrics are recorded as
plot failures. Engine-default PDF figures remain available alongside these
SVGs.

The generated `exploratory/configured_plots/` directory contains an SVG,
exact plotting-data TSV, and JSON settings recording coordinate/metadata
sources and axis limits. `plot_status.tsv` records every requested QC or
analysis plot, including invalid scopes, metrics, scales, metadata columns
and duplicate `plot_id` values. All rows sharing a duplicate ID fail;
requests are never silently discarded. A valid plot whose values or
coordinates are unavailable is `blocked`, with a reason. Analysis coordinates
may be unavailable while an analysis is pending, during validation-only runs,
or when normalization cannot run. `analysis_plot_status.tsv` retains the
analysis-only view of these outcomes.

## Explore before making comparisons

Set `analysis_mode=exploratory` for normalization and exploratory figures when
the experimental design or biological independence has not been confirmed.
Leave that analysis's contrast rows empty. Its output is descriptive and does
not produce DEG claims. Metadata supplied for display remains a declaration;
labels copied from an old sample table are not proof of treatment assignment.

```bash
scripts/symbiomics post validate -p my_study
scripts/symbiomics post run -p my_study
```

For inference, set `analysis_mode=inferential` (the default), declare the
experimental variables and confirmed biological units, and fill contrast
rows. Inspect the validation outputs before interpreting comparisons.

| Design | Configuration |
| --- | --- |
| Independent groups | `design_type=group`, `factor_cols=treatment`, explicit treatment reference |
| Additive effects | `design_type=additive`, declared factors/numeric predictors and optional batch/study nuisance columns |
| Paired observations | `design_type=paired`, `subject_col=subject_id` and within-subject factors |
| Factorial interaction | `design_type=factorial`, factors plus `interaction_terms=genotype:treatment` |
| Linear time/dose | `design_type=trend`, numeric predictor and declared unit; slope or specified change interval |

List columns with `|`. Declare all categorical reference levels, for example
`reference_levels=treatment=control;batch=run1`, and numerical units, for
example `numeric_units=time_h=hour`. Structured values containing reserved
`%`, `|`, `;`, `=`, `:` characters use UTF-8 percent encoding; ordinary
`samples.tsv` labels are retained unchanged. Table cells are never evaluated
as Python, shell or R code.

For a treatment versus control comparison, `contrasts.tsv` contains two rows:

```tsv
analysis_id	contrast_id	basis	term	weight	test_id
main	treated_vs_control	cell	treatment=treated	1	
main	treated_vs_control	cell	treatment=control	-1	
```

Cell terms specify every explanatory factor and numeric evaluation point.
Cell weights sum to zero. Coefficient contrasts use the public coefficient
IDs written by validation, such as `numeric:time_h`; a weight of `24` then
tests the modeled change over 24 declared time units. Do not copy internal R
coefficient names or rely on sample order.

The engine checks full model rank, residual degrees of freedom, estimability,
replication and repeated biological units. It reports confounding instead of
silently dropping batch or subject terms. Mixed models, general repeated
measures, nonlinear time models and single-cell inference are outside the
fixed-effect bulk-count interface.

Analysis-specific filtering and TMM precede edgeR quasi-likelihood analysis.
`qlf` with `lfc_cut` uses that effect-size cutoff to label/report results;
`effect_test=treat` tests an effect-size threshold directly. `effect_test=deseq2` runs DESeq2 (Wald test, TMM-derived size factors, ashr-shrunk log2FC, BH over all tested genes) on the same filtered genes, design matrix and contrast vectors; `logFC` then holds the shrunk effect, `log2FC_mle` the unshrunk one, and `F` is the squared Wald statistic. Omnibus (`test_id`) rows require `qlf`. BH correction is
within each tested analysis/contrast (or omnibus test), not a global FDR over
all comparisons. Batch-adjusted exploratory views, when requested and valid,
are not fed into DEG testing.

Within-target relative expression and the proportion of reads assigned to a
target answer different questions. This distinction matters particularly when
mainly host reads were aligned to a fungal reference. Sparse target signal can
support useful QC while leaving a requested DEG analysis mathematically
unavailable; the run preserves the QC and records the reason.

## Merge compatible execution projects

Merging is optional; `post run` uses the same importer internally.

```bash
scripts/symbiomics merge \
  --inputs post/my_study/inputs.tsv \
  --samples post/my_study/samples.tsv \
  --out /data/merged/my_study
```

Different assays remain separate. Within an assay, incompatible references,
annotations, counting policies/units, duplicate IDs or different gene sets
are errors. Gene rows are matched by ID; absent genes are not silently
zero-filled. Provenance records distinguish declared identifiers from verified
information and retain hashes of the actual input files.

Technical lane/run aggregation is explicitly requested with
`technical_merge_id`, `technical_partition_id`, `technical_disjoint=true` and
consistent library, biological-unit and experiment metadata. Its numerator
and denominator counts are summed before computing rates. The original QC
and aggregated QC are both retained. Different tissues/timepoints from one
subject must not be combined as technical partitions.

## Optional functional annotations

Supply `gene_go` tables with `gene_id, go_id` columns or `gene_pathway` tables
with `gene_id, pathway_id, pathway_name` columns, using tab delimiters.
`gene_labels` tables use `gene_id, label`. IDs remain case-sensitive; explicit
ID mappings are required when annotation IDs differ from the tested genes.

Request `enrichment=go` with `go_method=fisher_ora`, or
`enrichment=pathway` with `pathway_methods=camera`, `fgsea` or `camera|fgsea`.
Use `go|pathway` to request both. Each method reports its mapped background,
selection or ranking definition, and correction scope. Distinct pathway IDs
are retained even when their names match. Organism-specific gene exclusions
and pathway-pruning rules are not default settings.

## Runtime configuration

The post-analysis environment is described in `envs/post.yml`. Statistical
analysis requires Rscript and the engine's R/Bioconductor packages; QC-only
uses Python. The runtime selects an explicit `rscript` path first, then PATH,
then the `symbiomics-post` conda environment. Package versions and session
information from the actual runtime are recorded per run. An environment
specification is not a fully resolved package lock.

`config.yaml` is optional and uses safe YAML parsing (PyYAML). Supported keys:

```yaml
schema_version: 1
out_dir: output
threads: 1
strict_provenance: false
# rscript: /path/to/Rscript
timeout_seconds: 3600
```

`out_dir` is relative to the project unless absolute. Unknown settings are
rejected. Strict provenance requires reference and annotation hashes; a label
alone does not verify the corresponding reference bytes. Salmon inference
requires original quantification and transcript-to-gene mapping for tximport;
a preaggregated/rounded Salmon gene matrix or TPM matrix is not a substitute
for its effective-length-aware import.

## Outputs and partial failure

Each validate/run attempt creates a new directory under
`output/runs/<run_id>/`; reruns do not overwrite previous results. A lock
prevents concurrent writers. `latest_run.json` points to the latest completed
attempt and `latest_success.json` only changes for successful completion.
Resume/reuse of a previous run is not implemented.

| Path inside a run | Contents |
| --- | --- |
| `00_quant_qc/sample_metrics.tsv` | All observations, raw metrics, numerator/denominator, units, scope, availability reasons |
| `00_quant_qc/metric_definitions.tsv` | Metric definitions |
| `00_quant_qc/group_summary.tsv` | Descriptive summaries, valid/missing/zero counts, units and grouping metadata |
| `00_quant_qc/unassigned_long.tsv` | Original tool categories, counts and composition validity |
| `00_quant_qc/plots/` | Local SVG charts and plotting-data TSV files |
| `resolved/` | Resolved input/configuration tables and ID mappings |
| `analyses/<analysis_id>/design/` | Selected samples, coefficient/model/contrast preview and validation outputs |
| `analyses/<analysis_id>/normalized/` | Filtering, normalization and expression objects/tables |
| `analyses/<analysis_id>/exploratory/` | PCA/MDS, sample distances, heatmap and supporting data |
| `analyses/<analysis_id>/contrasts/` | Contrast results, volcano plots and requested enrichment |
| `stage_status.tsv`, `manifest.json` | Stage outcomes, reasons and execution provenance |
| `plot_status.tsv` | Every requested QC/analysis plot, including invalid requests and unavailable data |
| `analysis_plot_status.tsv` | Analysis-only view of requested PCA/MDS plot outcomes |
| `report.html` | Offline QC-first report and links to available results |

Statuses are `ok`, `skipped`, `failed` or `blocked`. Requested failures produce
a nonzero exit code even if QC or other analyses succeeded. Empty libraries,
no genes after filtering, unidentifiable models and unavailable R dependencies
retain available quantification QC. A QC-only run does not claim statistical
testing occurred. Optional unrequested analysis modules are recorded as
skipped rather than presented as empty successful discoveries.

## Run the post-analysis checks

From the repository root, the Python importer, CLI and report checks run with:

```bash
python3 -m unittest discover -s tests -p 'test_post_*.py'
```

The R acceptance cases are opt-in and otherwise reported as skipped. To run
them against an installed post-analysis R environment:

```bash
SYMBIOMICS_TEST_R=1 SYMBIOMICS_RSCRIPT=/path/to/Rscript \
  python3 -m unittest discover -s tests -p 'test_post_*.py'
```

These cases exercise the actual statistical engine and its dependencies;
set the Rscript path to the environment you intend to use for analysis.


## Multi-organism read composition: `post composition`

When reads were aligned separately to several genomes (host plus symbiont or
pathogen) and counted with HTSeq, `symbiomics post composition` reports how the
input read pairs split between genomes.

```bash
symbiomics post composition --libraries libraries.tsv --out composition \
    --primary host --signal-organism symbiont --signal-ppm 1000 --min-aligned-pairs 5e6
```

`libraries.tsv` has one row per library and genome: `sample_id`, `organism`,
`file` (raw two-column HTSeq output with `__*` special counters; relative paths
are resolved against the table). The `--primary` organism's counters define the
input pair total; each organism adds aligned pairs (total minus `__not_aligned`),
aligned percent, assigned pairs and assigned pairs per million input pairs.
Outputs: `composition.tsv`, `composition.svg` (stacked bars, remainder grey).

`--signal-organism` adds `signal_group` (`high` at or above `--signal-ppm`, else
`low`); `--min-aligned-pairs` adds `qc_pass` on the primary organism. Signal
groups describe read abundance only. They are not an infection rate or dose and
must not be correlated with expression as one. Genomes are aligned independently,
so the per-organism percentages can overlap (multi-mapping across genomes).
