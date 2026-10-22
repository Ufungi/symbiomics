# Usage

## Run it

```bash
scripts/symbiomics preflight                       # check the environment first

scripts/symbiomics run . -profile singularity,local64 \
    --input samplesheet.csv \
    --genome genome.fasta \
    --project plant
```

`scripts/symbiomics` is the supported entry point. Calling `nextflow` directly
works only if your `java` is already 17–21; the launcher exists because many
hosts default to Java 11, and it points `JAVA_CMD` at the JDK inside the
Nextflow conda env.

### Plan first

```bash
scripts/symbiomics run . -entry strategy --genome genome.fasta --taxon plant
```

Measures the genome, resolves every decision, prints `strategy.yml`, runs
nothing else. Do this before launching anything that will take days.

### Resume

```bash
scripts/symbiomics run . -profile singularity,local64 ... -resume
```

### Subsets

```bash
--steps align,assemble          # only these stages
--skip_steps quantify           # everything except this
--sample MySE01               # one row of the samplesheet
--test_mode --test_reads 200000 # same code paths, less data
```

`--test_mode` never changes which code path runs — same masker, same aligner,
same decisions — only the data volume. That is the point of it.

## Existing counts: quantification QC and post-analysis

```bash
scripts/symbiomics post init -p my_study
# Fill post/my_study/inputs.tsv and samples.tsv.
scripts/symbiomics post run -p my_study --qc-only
# After completing and confirming the experimental design tables:
scripts/symbiomics post validate -p my_study
scripts/symbiomics post run -p my_study
```

The `post` command family is standalone and does not launch Nextflow. It
preserves zero/low-signal samples in QC and records separate outcomes for
requested analyses. Use `analysis_mode=exploratory` for PCA before a design is
confirmed. See [post-analysis](post_analysis.md) for table contracts,
`merge`, explicit project directories, runtime configuration and limitations.

## The samplesheet

CSV by default; the template (`assets/samplesheet.example.csv`) opens directly
in Excel — a `#` comment block at the top explains each column, followed by one
example row per layout (SE / PE / pre-aligned BAM). Only `sample_id` plus one of
`fastq_1` / `bam` is required.

```csv
sample_id,fastq_1,fastq_2,strandedness,bam
MySE01,/raw/MySE01.fastq.gz,-,auto,-
MyPE01,/raw/MyPE01_R1.fastq.gz,/raw/MyPE01_R2.fastq.gz,auto,-
MyBAM,-,-,reverse,/aligned/P.sorted.bam
```

| column | values | default |
|---|---|---|
| `sample_id` | `[A-Za-z0-9][A-Za-z0-9._-]{0,63}`, unique | required |
| `layout` | `SE`, `PE`, `interleaved`, `auto`, `-` | `auto` |
| `fastq_1` / `fastq_2` | path, or `;`-separated list of technical replicates | required / `-` |
| `strandedness` | `auto`, `unstranded`, `forward`, `reverse` | `--strandedness` |
| `read_type` | `short`, `isoseq`, `ont_cdna`, `ont_drna` | `short` |
| `bam` | path to a coordinate-sorted, indexed BAM | `-` |
| `use_for` | subset of `braker,stringtie,count`, or `all` | `all` |

Lines starting with `#` and blank lines are ignored. Relative paths resolve
against `--raw_dir`. Header aliases are accepted (`sample`, `R1`, `read1`, `fq1`…).

`use_for` is how you say "these samples feed StringTie and counting, but only
some of them are gene-prediction evidence" without maintaining two sheets.

A resolved copy is always written to
`project/<project>/output/00_pipeline_info/samplesheet.resolved.tsv`. That file is the
record of what actually ran and can be fed straight back in as `--input`.

### SE/PE detection

An explicit `layout` always wins. Otherwise:

1. `fastq_2` set → PE.
2. Otherwise the mate is looked for on disk: `_R1_001`→`_R2_001`, `_R1`→`_R2`,
   `.R1`→`.R2`, `_1`→`_2`, `.1`→`.2`, `_fwd`→`_rev`, `_f`→`_r`, `-1`→`-2`.
   A candidate is only accepted if it exists, is non-empty, and its first four
   read IDs match the R1 file's.
3. Nothing found → SE.

`interleaved` is never inferred, only declared.

A sheet that contradicts itself fails rather than being corrected. Declaring
`SE` while `fastq_2` holds a file is an error — silently discarding half a
library is not a recoverable mistake.

## Strandedness

`auto` (the default) infers it per sample with no annotation, by reading splice
junction motifs out of the genome: `GT..AG` means the transcript is on the plus
strand, `CT..AC` means minus, and the fraction that agrees with the read's own
strand is the answer.

Salmon- and RSeQC-based inference need a transcriptome or a BED12 gene model,
which do not exist yet at this point in the pipeline — producing them is why the
pipeline is running. Junction motifs need only the genome.

| f | call |
|---|---|
| < 20,000 informative reads | `undetermined` |
| ≥ 0.85 | `forward` |
| ≤ 0.15 | `reverse` |
| 0.40–0.60 | `unstranded` |
| anything else | `ambiguous` |

`ambiguous` is deliberately separate from `unstranded`. A value of 0.72 is not a
weakly stranded library, it is a signal of mixed prep, gDNA carry-over, or
degraded RNA, and calling it `forward` would hide that.

A declared value in the samplesheet always wins, but inference still runs and a
disagreement is reported. `--strict_strandedness` turns disagreements and
undetermined calls into errors.

If samples disagree with each other (some forward, some reverse), the pipeline
warns loudly: mixing them in one differential-expression analysis is not valid.

## Output

`project` (in `config.yaml`) or `--project NAME` (CLI) decides the output
directory: everything lands in `project/<project>/output/`.

```
project/<project>/output/
├── 00_pipeline_info/   strategy.yml, samplesheet.resolved.tsv, versions.yml, traces
├── 10_genome/          .fai, genome stats
├── 30_reads/           fastqc, trimming reports
├── 35_strandedness/    per-sample calls with the evidence behind them
├── 40_align/           *.sorted.cram (+ .crai), summaries, samtools stats
├── 45_assemble/        per-sample and merged StringTie GTF
├── 95_quantify/        counts_featurecounts.tsv, counts_htseq.tsv, per-sample QC
└── 99_report/          multiqc_report.html
```

### Publishing large alignment sets

Alignments are stored as CRAM (reference-compressed) and published as
`*.cram` + `*.crai` when `--save_bam` is on (the default). `--publish_mode copy`
is the default. `link` (hardlink) avoids duplicating the CRAMs but only works
when the Nextflow work directory and the project output
(`project/<project>/output/`) are on the **same filesystem**; the launcher puts work on
`/data/db` by default, so set it explicitly first:

```bash
NXF_WORK=/path/to/work scripts/symbiomics run . ... --publish_mode link
```

The pipeline checks this at startup and refuses rather than failing after the
first large CRAM is written.

## Profiles

| profile | effect |
|---|---|
| `singularity` | containers (recommended) |
| `docker` | containers via Docker |
| `conda` | conda envs; see the gap list in `conf/conda.config` |
| `gpu` | `--nv` passthrough, GPU tasks serialised |
| `local64` | 64-core / 480 GB host with no scheduler |
| `test` | 75 kb smoke test |
| `test_fungus` | *T. matsutake*, 161 Mb |

Combine them: `-profile test,singularity`.

## Reusing an index

```bash
--hisat2_index /data/.../P.densiflora_v1.0_HA.genome
```

A **prefix**, in hisat2's sense — the part before `.1.ht2`, not a single file.
The pipeline verifies the index exists before starting and reports how many
files it found. On a 22 Gb genome this saves days.
