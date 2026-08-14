# Usage

## Run it

```bash
scripts/symbiomics preflight                       # check the environment first

scripts/symbiomics run . -profile singularity,local64 \
    --input samplesheet.tsv \
    --genome genome.fasta \
    --outdir results
```

`scripts/symbiomics` is the supported entry point. Calling `nextflow` directly
works only if your `java` is already 17–24; the launcher exists because many
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
--sample SRR8356879             # one row of the samplesheet
--test_mode --test_reads 200000 # same code paths, less data
```

`--test_mode` never changes which code path runs — same masker, same aligner,
same decisions — only the data volume. That is the point of it.

## The samplesheet

TSV by default; `.csv` is detected by extension. Only `sample_id` plus one of
`fastq_1` / `bam` is required.

```tsv
sample_id	layout	fastq_1	fastq_2	strandedness	condition	replicate	read_type	bam	genome_id	use_for
SRR8356879	SE	SRR8356879.fastq.gz	-	auto	fruitbody	1	short	-	Tmat	all
MyPE01	-	MyPE01_R1.fq.gz	MyPE01_R2.fq.gz	auto	mycelium	1	short	-	Tmat	all
MyPE02	-	A_R1.fq.gz,B_R1.fq.gz	A_R2.fq.gz,B_R2.fq.gz	forward	mycelium	2	short	-	Tmat	all
Prealign	-	-	-	reverse	needle	1	short	/x/P.sorted.bam	Pinde	stringtie,count
```

| column | values | default |
|---|---|---|
| `sample_id` | `[A-Za-z][A-Za-z0-9._-]{0,63}`, unique | required |
| `layout` | `SE`, `PE`, `interleaved`, `auto`, `-` | `auto` |
| `fastq_1` / `fastq_2` | path, or comma list of technical replicates (`;` in CSV) | required / `-` |
| `strandedness` | `auto`, `unstranded`, `forward`, `reverse` | `--strandedness` |
| `read_type` | `short`, `isoseq`, `ont_cdna`, `ont_drna` | `short` |
| `bam` | path to a coordinate-sorted, indexed BAM | `-` |
| `use_for` | subset of `braker,stringtie,count`, or `all` | `all` |

Lines starting with `#` and blank lines are ignored. Relative paths resolve
against `--raw_dir`. Header aliases are accepted (`sample`, `R1`, `read1`, `fq1`…).

`use_for` is how you say "these samples feed StringTie and counting, but only
some of them are gene-prediction evidence" without maintaining two sheets.

A resolved copy is always written to
`results/00_pipeline_info/samplesheet.resolved.tsv`. That file is the record of
what actually ran and can be fed straight back in as `--input`.

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

```
results/
├── 00_pipeline_info/   strategy.yml, samplesheet.resolved.tsv, versions.yml, traces
├── 10_genome/          .fai, genome stats
├── 30_reads/           fastqc, trimming reports
├── 35_strandedness/    per-sample calls with the evidence behind them
├── 40_align/           *.sorted.bam (+ .csi/.bai), summaries, samtools stats
├── 45_assemble/        per-sample and merged StringTie GTF
├── 95_quantify/        counts_featurecounts.tsv, counts_htseq.tsv, per-sample QC
└── 99_report/          multiqc_report.html
```

### Publishing large BAM sets

`--publish_mode copy` is the default. `link` (hardlink) avoids duplicating the
BAMs but only works when the Nextflow work directory and `--outdir` are on the
**same filesystem**; the launcher puts work on `/data/db` by default, so set it
explicitly first:

```bash
NXF_WORK=/path/to/results/work scripts/symbiomics run . ... --publish_mode link
```

The pipeline checks this at startup and refuses rather than failing after the
first large BAM is written.

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
