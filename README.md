# eukannot

A decision-driven Nextflow pipeline for eukaryote genome annotation and mRNA-seq
quantification. Built to work on a 160 Mb fungal genome and a 22 Gb conifer
genome with the same command.

```
scripts/eukannot run . -profile singularity,local64 \
    --input samplesheet.tsv --genome genome.fasta --outdir results
```

## What makes it different

Most annotation pipelines make you pick the tools. This one **decides, records
why, and degrades gracefully**.

After reading the genome, a policy engine writes `strategy.yml`:

```yaml
decisions:
  masker:
    choice: red
    reason: "size_class=huge, taxon=plant; de-novo modelling does not finish at
             this scale (measured: 52 h to an empty library on a 21.7 Gb conifer)"
    alternatives: [repeatmodeler_subsample]
    escape: ["--repeat_lib <fasta>", "--premasked", "--masker <name>"]
    blocked_at_this_size: [earlgrey, edta]
  bam_index:
    choice: csi
    reason: "longest contig 2,007,914,973 bp exceeds the BAI ceiling 536,870,912 bp"
budget:
  braker_sharded:
    max_time: 240.h
    on_exceed: continue_without
    fallback: "drop braker_sharded; finish with the remaining tracks and record
               the omission in the report"
```

Three properties follow from that:

1. **Nothing is a black box.** Every choice carries its reason, the alternatives
   it rejected, and the escape hatches that override it.
2. **You can review the plan before committing.** `--dry_run_strategy` resolves
   everything and stops. A 22 Gb genome is a multi-week run; the plan should be
   read first, not reconstructed from the log afterwards.
3. **Partial failure is not total failure.** Tracks carry time budgets. One that
   blows its budget is dropped, the pipeline finishes with what is left, and the
   omission is recorded in the report rather than silently ignored.

Requesting something that cannot work fails immediately with the arithmetic:

```
ERROR ~ [align] bam_index=bai requested but the longest contig is 2,007,914,973 bp,
        above the BAI limit of 536,870,912.
        fix: --bam_index csi
```

## Status

**v0.1.0 — the mRNA-seq arm.** Working end to end:

- samplesheet parsing with SE/PE auto-detection
- annotation-free strandedness inference
- HISAT2 alignment (large-index and mmap aware), coordinate sort, CSI/BAI index
- StringTie assembly and merge
- featureCounts / HTSeq counting, merged matrices, MultiQC
- the policy engine and `strategy.yml`

Structural annotation (repeats, BRAKER, Helixer, consensus), functional
annotation, and the sharded huge-genome path land in later milestones. See
`docs/roadmap.md`.

## Install

```bash
git clone https://github.com/Ufungi/eukannot.git
cd eukannot
scripts/eukannot preflight        # checks java, nextflow, containers, GPU, disk
```

`scripts/eukannot` is the supported entry point. It exists because Nextflow needs
Java 17–24 and many hosts default to Java 11; the launcher points `JAVA_CMD` at
the JDK bundled with the Nextflow conda env.

## Input

One TSV (or CSV). Only `sample_id` and one of `fastq_1` / `bam` are required.

```tsv
sample_id	layout	fastq_1	fastq_2	strandedness	condition	replicate	read_type	bam	genome_id	use_for
SRR8356879	SE	SRR8356879.fastq.gz	-	auto	fruitbody	1	short	-	Tmat	all
MyPE01	-	MyPE01_R1.fq.gz	MyPE01_R2.fq.gz	auto	mycelium	1	short	-	Tmat	all
Prealign	-	-	-	reverse	needle	1	short	/x/P.sorted.bam	Pinde	stringtie,count
```

`layout` may be left as `-`: the mate is inferred from the filesystem
(`_R1`→`_R2`, `_1`→`_2`, …) and the pairing is confirmed against the read IDs.
An explicit `layout` always wins, and a sheet that contradicts itself fails
loudly rather than being quietly corrected — declaring `SE` while a mate file
sits in `fastq_2` is an error, not a silent discard.

See `assets/samplesheet.example.tsv`, `assets/genomes.example.tsv` (multi-genome
batches) and `assets/reference_proteomes.example.tsv`.

## Docs

| | |
|---|---|
| `docs/usage.md` | running it, parameters, output layout |
| `docs/decisions.md` | how the policy engine chooses, and how to override it |
| `docs/huge_genomes.md` | the 22 Gb conifer path: masking, sharding, CSI, index reuse |
| `docs/clade_notes.md` | which tools do not support which clades, with citations |
| `docs/troubleshooting.md` | the errors you will actually hit |

## License

MIT.
