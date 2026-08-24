# Multi-genome symbiont mapping (`-entry symbiont_mapping`)

## Why

Every other entry point in this pipeline aligns reads against exactly one
genome, decided before the run starts. That is the right default for a
clean single-organism sample, but it does not fit the actual research
subject this pipeline is named after: a mycorrhizal root or fruiting-body
sample carries reads from *both* the host (*Pinus densiflora*) and the
ectomycorrhizal symbiont (*Tricholoma matsutake*), or a sample's symbiont
identity is one of several candidate species rather than known up front.
Committing such a sample to one reference before alignment either mixes the
other organism's reads in as noise or throws them away as "unmapped" with no
record of where they might actually belong.

`-entry symbiont_mapping` aligns the same trimmed reads against every genome
in `--mapping_genomes` and reports, per sample, which genome the reads
actually match -- and how confidently.

## Usage

```bash
scripts/symbiomics run . -entry symbiont_mapping -profile singularity,local64 \
    --input samplesheet.tsv \
    --mapping_genomes mapping_genomes.tsv \
    --outdir results_mapping
```

`--mapping_genomes` is a TSV, one row per candidate genome (see
`assets/mapping_genomes.example.tsv`):

| column | required | meaning |
|---|---|---|
| `genome_id` | yes | unique short id, used in output filenames and the comparison table |
| `fasta` | yes | path to the genome FASTA |
| `taxon` | no | informational only (`auto`\|`fungi`\|`plant`\|`animal`\|`other`) |

At least 2 genomes are required -- a single genome has nothing to be
compared against; use the default entry point for that case instead.

## What runs

Deliberately simpler than the main align stage: no strandedness probe pass.
`--rna-strandness` mainly steers spliced-alignment sensitivity and the XS tag
StringTie needs, not the overall alignment rate this comparison is built on,
so a sample's declared (or default `auto` -> unstranded-equivalent) value is
passed straight through HISAT2 rather than spending an extra probe alignment
per sample to infer it. If a sample declares `forward`/`reverse` in the
samplesheet that is honoured as-is.

1. `fastp` trims each sample once (not once per genome).
2. A HISAT2 index is built once per candidate genome.
3. Every sample is aligned against every genome: N samples x G genomes =
   N*G alignment tasks. This is the deliberate cost of the comparison --
   there is no shortcut that scores a sample against a genome without
   aligning it.
4. `samtools sort`/`index`/`stats`/`flagstat`/`idxstats` per (sample, genome)
   pair, same as the main align stage.
5. `bin/summarize_multi_genome_mapping.py` collects every HISAT2
   `--new-summary` overall alignment rate into one matrix and picks, per
   sample, the genome with the highest rate.

## Output (`60_symbiont_mapping/`)

- `mapping_matrix.tsv` -- sample_id x genome_id, overall HISAT2 alignment
  rate (%) for every combination.
- `best_genome.tsv` -- per sample: `best_genome`, `best_rate`,
  `runner_up_genome`, `runner_up_rate`, `margin`, `ambiguous`.
- `mapping_summary.json` -- the same information as structured JSON, plus
  any summary files that failed to parse.
- `index/`, `log/`, `stats/` and the sorted/indexed BAMs themselves (subject
  to `--save_bam`), one set per (sample, genome) pair.

## Reading `ambiguous`

A sample is flagged `ambiguous` when its top two genomes' alignment rates
are within `--mapping_ambiguous_margin` percentage points (default `5.0`).
This is a coverage argmax, not a taxonomic identification: `ambiguous` on
every sample of a
project usually means the candidate genomes are too closely related for
short-read HISAT2 alignment rate alone to separate them (near-identical
paralogs, an unremoved shared organelle contig), and a uniformly low rate on
every candidate means the reads do not belong to any genome offered, not
that the call is a coin flip between them.

## Verification status

`bin/summarize_multi_genome_mapping.py` is unit-tested
(`tests/test_summarize_multi_genome_mapping.py`) against synthetic HISAT2
summary files. The Nextflow side (`modules/local/multi_genome_map.nf`,
`subworkflows/local/multi_genome_mapping`, the `symbiont_mapping` entry
point in `main.nf`) was written and hand-reviewed against this pipeline's
existing DSL2 conventions, but has **not** been run: this was built in a
session with no Nextflow runtime available, so neither the added
`-entry symbiont_mapping -stub-run` check in `tests/run_tests.sh` nor a real
alignment has actually been executed. Run
`scripts/symbiomics run . -entry symbiont_mapping -profile test_mapping -stub-run`
(DAG wiring) and then a real run with real candidate genomes before relying
on this entry point.

## What this deliberately is not

This is read-count triage by alignment rate, not deconvolution: it does not
attempt to resolve reads that align equally well to two genomes down to the
individual read, and it makes no attempt at species identification beyond
"best match among the genomes you supplied." If the intended use is producing
a per-sample host and symbiont BAM to feed into `QUANTIFY` downstream, take
`best_genome.tsv`'s call, then rerun the default entry point per organism
with `--sample` scoped to the samples assigned to it -- this entry point's
job ends at the comparison, deliberately, rather than silently forking into
a second full pipeline run.
