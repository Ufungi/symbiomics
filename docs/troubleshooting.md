# Troubleshooting

Errors listed roughly in the order you are likely to hit them.

## `Cannot find Java or it's a wrong version`

Nextflow needs Java 17–24. Many hosts default to Java 11.

```
NOTE: Nextflow is trying to use the Java VM defined by the following environment variables...
```

**Fix:** use `scripts/eukannot`, which points `JAVA_CMD` at the JDK bundled with
the Nextflow conda env. If you must call `nextflow` directly:

```bash
conda activate nextflow
export JAVA_HOME="$CONDA_PREFIX/lib/jvm"
export JAVA_CMD="$JAVA_HOME/bin/java"
```

`scripts/eukannot preflight` checks this first and prints the exact remedy.

## `Failed to publish file ... [link]`

```
Failed to publish file: /data/db/.../mini.fasta.fai; to: /data/snyoo/.../10_genome/mini.fasta.fai [link]
```

Hardlinks cannot cross filesystems, and the launcher puts the work directory on
`/data/db` while `--outdir` is usually elsewhere.

**Fix:** either keep the default `--publish_mode copy`, or put the work dir on
the same filesystem as the output:

```bash
NXF_WORK=/path/to/results/work scripts/eukannot run . ... --publish_mode link
```

The pipeline now checks this at startup and refuses up front rather than failing
after the first large BAM has been written.

## `Failed to pull singularity image ... the requested image was not found`

A container tag that does not exist. Biocontainer tags look predictable but are
not always published.

**Fix:**

```bash
scripts/check_containers.sh
```

This HEAD-checks every image in `conf/containers.config` and names the ones that
do not resolve. Run it in CI weekly and before tagging a release.

## `mkfifo(/tmp/NN.inpipe1) failed`

HISAT2 creates FIFOs to decompress gzipped reads, and `/tmp` inside a container
is frequently read-only.

**Fix:** already handled — `HISAT2_ALIGN` exports `TMPDIR="$PWD"`. If you see
this in a module you have added, do the same.

## `exit status 127` in a piped process

A process that pipes two tools together needs both in one image. The plain
`hisat2` biocontainer has no `samtools`, which is why `HISAT2_ALIGN` uses a
mulled image carrying both.

**Fix:** point the process at a mulled container, or add both tools to its conda
env. Do not "solve" it by writing an intermediate SAM to disk.

## `Cannot load from object array because "this.keys" is null`

A `groovy.json.JsonSlurper` `LazyMap` being read from several dataflow actors at
once. The map defers building its backing store and the deferred build races.

**Fix:** parse with `JsonSlurperClassic`, which returns plain `LinkedHashMap`s.
`GENOME_PREP` does this for `strategy.json`.

## A process never starts and the run hangs, then aborts

Check the "still active" block at the end of `.nextflow.log`:

```
[process] EUKANNOT:RNASEQ_ALIGN:HISAT2_BUILD
  status=ACTIVE
  port 1: (value) OPEN  ; channel: large_index
```

`OPEN` means that input never received a value. Usually one of:

- a **queue channel consumed twice** — the second consumer starves. Add
  `.first()` to turn it into a value channel when several processes read it.
- an upstream `map` closure threw, so nothing was ever emitted.

## `[samplesheet] ...` errors

These are deliberate and each names the line, what was tried, and the fix.

| message | meaning |
|---|---|
| `layout=SE but fastq_2=... exists` | the sheet contradicts itself; half the library would be discarded |
| `layout=PE but fastq_2 is '-'` | no mate given |
| `inferred mate ... but read IDs disagree` | the filesystem mate is not actually the mate |
| `fastq not found. Tried: ...` | shows the path after `--raw_dir` resolution |
| `sample_id ... is duplicated` | use the `replicate` column instead |
| `bam is name-sorted (SO:queryname)` | coordinate-sorted BAM is required |
| `bam has no index` | run `samtools index -c` |

## `[repeats] --premasked was given but ... is 0.00% lowercase`

Soft-masking is case-based; a genome with no lowercase carries no mask. Either
drop `--premasked` and let the pipeline mask it, or supply a genuinely
soft-masked genome (`RepeatMasker -xsmall`).

## `[repeats] masker='edta' cannot complete on a huge genome`

De-novo repeat modelling does not finish at conifer scale — see
`docs/huge_genomes.md` for the measurements. Use `--masker red`,
`--masker repeatmodeler_subsample`, `--repeat_lib`, or `--premasked`.

## `[align] bam_index=bai requested but the longest contig is ...`

BAI cannot address a contig above 512 Mb. Use `--bam_index csi`. There is no
workaround; the format does not support it.

## Strandedness is `undetermined` or `ambiguous`

- `undetermined`: fewer than `--strandedness_min_reads` (20,000) informative
  spliced reads. Common on very small test data or an intron-poor genome.
- `ambiguous`: the agreeing fraction fell between the stranded and unstranded
  bands. Treat this as a signal — mixed library prep, gDNA carry-over, or
  degraded RNA — not as noise.

Both fall back to `unstranded` with a warning, or fail under
`--strict_strandedness`. To skip inference entirely, set `--strandedness`
globally or per row in the samplesheet.

## Warnings that are safe to ignore

```
WARN: There's no process matching config selector: TRIMGALORE
WARN: There's no process matching config selector: AGAT_.*
```

`conf/modules.config` and `conf/containers.config` carry entries for stages that
are not yet in the v0.1 DAG. They become live in later milestones.

## Getting more detail

```bash
scripts/eukannot run . ... -with-trace -with-report -with-timeline
less .nextflow.log
cd <work dir printed in the error>; cat .command.sh .command.err
```
