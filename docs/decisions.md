# How the pipeline decides

The centre of eukannot is not a tool list, it is a policy engine. After the
genome is measured, `DECIDE_STRATEGY` writes `strategy.yml` and every stage
downstream reads it instead of re-deriving the same choice in five places.

## Inputs to the decision

| Input | Where it comes from |
|---|---|
| genome size, contig count, longest contig, N%, softmask % | `bin/genome_stats.py`, one streaming pass |
| taxon / clade | `--taxon`, `--clade`, or `--species_taxid` |
| available evidence | the samplesheet (RNA / long reads) and `--proteins` |
| resources | `--max_cpus`, `--max_memory`, GPU presence |
| user budget | `--time_budget`, `--priority` |
| explicit overrides | any parameter the user set |

## Output

```yaml
decisions:
  <name>:
    choice: <what will run>
    reason: <why, in terms of the measured inputs>
    source: policy | user
    alternatives: [...]          # what was rejected but is still reachable
    escape: [...]                # the flags that override this
    blocked_at_this_size: [...]  # what cannot work here at all
budget:
  <track>:
    max_time: 240.h
    on_exceed: continue_without
    fallback: <what the pipeline does instead>
notes: [...]                     # anything the user should know
```

## Review before you commit

```bash
scripts/eukannot run . -entry strategy --genome big.fasta --taxon plant
```

Runs the measurement and the decision, prints the plan, and stops. For a
multi-week genome this is the difference between reviewing the plan and
reconstructing it from a log afterwards.

## The rules

### Genome size class

| class | size | example |
|---|---|---|
| tiny | ≤ 10 Mb | organelles, smoke tests |
| small | ≤ 500 Mb | *Tricholoma matsutake*, 161 Mb |
| medium | ≤ 2 Gb | rice, tomato |
| large | ≤ 5 Gb | wheat-scale |
| huge | > 5 Gb | *Pinus densiflora*, 21.74 Gb |

Thresholds are `--size_threshold_{tiny,small,medium,large}`.

### Repeat masking

| class | default | blocked |
|---|---|---|
| tiny | `red` | |
| small | `earlgrey` | |
| medium | `edta` (plant) / `earlgrey` | |
| large | `repeatmodeler_subsample` | `earlgrey`, `edta` |
| huge | `red` | `earlgrey`, `edta` |

The blocks are not conservatism. EarlGrey's own guidance is "weeks" at 25 Gb,
and a RepeatModeler run on the 21.7 Gb pine genome in this lab ran 52 hours
through round 4 and emitted a zero-byte library. Asking for a blocked masker is
an error that names the alternatives.

Escapes: `--repeat_lib <fasta>`, `--premasked`, `--masker <name>`.

`--premasked` is verified, not trusted: a genome declared premasked that is less
than 2% lowercase is rejected, because soft-masking is case-based and a genome
with no lowercase carries no mask. The check cannot prove the converse — some
assemblers emit lowercase of their own — so a high fraction is necessary but not
sufficient evidence that RepeatMasker ever ran.

### Aligner

STAR is refused above `--star_max_genome` (2 Gb default): its suffix array is
not buildable at conifer scale. HISAT2 gains `--large-index` above 4 Gb and
`--mm` once the index exceeds 8 GB, so concurrent aligners share one page-cache
copy instead of holding one resident copy each.

### BAM index

CSI whenever any contig exceeds 512 Mb, because BAI cannot address one at all —
the 2.008 Gb pine chromosome is the case that forces this. Below that limit both
are emitted, since some downstream tools still hardcode a `.bai` suffix.

### Evidence mode

| RNA | protein | IsoSeq | mode |
|---|---|---|---|
| ✓ | ✓ | ✓ | `dual` |
| ✓ | ✓ | | `etp` |
| ✓ | | | `et` |
| | ✓ | | `ep` |
| | | ✓ | `isoseq` |
| | | | `es` + warning |

`--evidence_mode ep` with RNA present is legal and means "ignore RNA for gene
prediction"; the BAMs still drive StringTie and counting.

### Predictor tracks

Up to `large`, BRAKER is primary. At `huge` the primary becomes Helixer with
BRAKER following in sharded form, because BRAKER has no internal checkpointing
and is not guaranteed to finish at 22 Gb, while Helixer completes in days on one
GPU. Both feed the consensus. With no GPU, this falls back to sharded BRAKER and
says so.

Consensus follows from the tracks: one track → none; BRAKER-family only →
TSEBRA; anything with non-BRAKER isoforms → Mikado.

### Clade exclusions

Some tools have no model for some clades. These are recorded rather than
rediscovered — see `docs/clade_notes.md`.

| tool | excluded | why |
|---|---|---|
| Tiberius | gymnosperms | ships Angiosperms (Mesangiospermae) only |
| EGAPx | gymnosperms, fungi | supports Magnoliopsida; fungi explicitly out of scope |

### Budgets and graceful degradation

Every heavy track carries a wall-clock budget. A track that exceeds it is
dropped, the pipeline finishes with the remaining tracks, and the omission is
written into `strategy.yml` and the final report. A pipeline that runs for two
weeks and then returns nothing because one optional track hung is a worse
outcome than one that returns a Helixer annotation and says BRAKER was dropped.

Override globally with `--time_budget 120.h`.

## Overriding

Every decision has an explicit parameter. Setting it flips `source` to `user`
and the reason records that. Overrides that cannot work still fail:

```
ERROR ~ [align] bam_index=bai requested but the longest contig is 2,007,914,973 bp,
        above the BAI limit of 536,870,912.
        fix: --bam_index csi
```

This is deliberate. A silent downgrade to something that "works" would produce
an index that cannot address half the genome.
