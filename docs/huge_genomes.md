# Huge genomes

Reference case: ***Pinus densiflora***, **21.74 Gb**, 2,006 contigs, longest
contig **2,007,914,973 bp**, **not softmasked**.

Almost every number in `conf/genome_size/huge.config` and every rule in the
`huge` branch of the policy engine exists because of something measured on this
genome. This page records what and why.

## What breaks at this scale

| | measurement | consequence |
|---|---|---|
| **De-novo repeat modelling** | a RepeatModeler run on this genome ran **52 h through round 4** and emitted a **0-byte `consensi.fa`**; only `round-1/consensi.fa` (1.19 MB) held anything. EarlGrey's own guidance is "weeks" at 25 Gb | EarlGrey and EDTA are **blocked**; default is Red |
| **STAR** | index build is not feasible at 21.7 Gb | STAR is **refused**; HISAT2 with `--large-index` |
| **BAI** | caps at 512 Mb per contig; the longest here is 2.008 Gb | **CSI only** |
| **HISAT2 index** | 38 GB on disk; N resident copies × 13 forks exceeds 480 GB | `--mm` forced, `maxForks 4` |
| **BRAKER** | no internal checkpointing; multi-week with no partial result | **sharded**, and Helixer runs first |
| **Existing GFF3** | `exon`/`CDS` carry `ID=` equal to the gene, no `Parent=`, 44,233 genes but only 4,637 `mRNA`, `CDS` and `cds` both present | must be repaired with AGAT before any use as a guide or comparison |

## Repeat masking

Three paths, in increasing cost:

### 1. Red (default)

```bash
--masker red
```

Machine-learning masker, no library, no de-novo modelling. Hours, not weeks.
This is also what BRAKER4 would do internally. Accept it unless you need TE
biology, not just a mask.

### 2. Model on a stratified subsample

```bash
--masker repeatmodeler_subsample --rm_sample_size 500000000
```

```
sample_contigs.py --total 500000000 --window 5000000 --min-contig 1000000 \
                  --exclude-n-frac 0.5 --seed 42 --out subsample.fa
```

Windows are spread **proportionally across all contigs**, not taken from the
longest ones. On this genome the single longest contig is 2 Gb and carries
nothing like the full TE diversity, so "take the biggest scaffolds" is the wrong
sampling strategy.

`-LTRStruct` is **off** on the subsample — it is the stage that consumed the 52
hours and produced nothing. Structural LTR discovery is re-added afterwards and
far more cheaply with LTR_retriever / LTRharvest on the same subsample.

The full genome is then masked with the subsample-derived library in 100 Mb
chunks that respect contig boundaries, 8 concurrent tasks × `-pa 8` = 64 cores.

### 3. Bring your own

```bash
--repeat_lib /path/to/families.fa       # skip modelling entirely
--premasked                             # genome is already softmasked
```

`--premasked` is checked, not trusted. This genome is 0.00% lowercase, so the
check rejects it:

```
ERROR ~ [repeats] --premasked was given but Pinde is 0.00% lowercase, i.e. not softmasked.
        Soft-masking is case-based; a genome with no lowercase carries no mask.
```

### Salvage

If a RepeatModeler run dies, `salvage_repeatmodeler.py` concatenates
`round-*/consensi.fa` and clusters with cd-hit-est. Be realistic about the
yield: in the observed failure only round 1 held anything, so this recovers a
partial library, not a good one. It is a fallback, not a plan.

## Alignment

```bash
--hisat2_index /data/.../P.densiflora_v1.0_HA.genome
```

Reuse a prebuilt index — the `.ht2l` set for this genome already exists and
building it again costs days. Pass the **prefix** (the part before `.1.ht2l`),
not a single file. The pipeline verifies it and reports how many index files it
found before starting.

`--mm` is forced above an 8 GB index so concurrent aligners share one page-cache
copy. Without it, four aligners against a 38 GB index is 152 GB of resident
memory instead of 38.

## BRAKER sharding

```
TRAIN    ~300 Mb of gene-dense contigs (ranked by StringTie exon coverage,
         excluding the top 1% by repeat density) -> AUGUSTUS species params
         + GeneMark model, using all evidence restricted to that subset

PREDICT  10 Mb windows, 500 kb overlap -> ~2,170 shards
         500 kb because conifer introns routinely exceed 100 kb
         augustus --predictionStart/--predictionEnd + gmes_petap --predict_only

MERGE    keep a gene iff its MIDPOINT falls in the window's core (window minus
         overlap); genes spanning more than the overlap go to boundary_genes.tsv
         and are re-predicted on merged windows
```

**Sharding is about checkpointing, not wall-clock.** It converts one
un-resumable multi-week process into 2,170 independently resumable tasks. That
is the difference between a crash costing an hour and costing a fortnight.

## Predictor order

At `huge` the policy engine runs **Helixer first, BRAKER second**:

- Helixer (`land_plant`) finishes in days on one RTX 3090 and yields a usable
  annotation while BRAKER is still running.
- BRAKER follows in sharded form and both feed the consensus.
- If BRAKER exceeds its budget it is dropped, the run completes on Helixer
  alone, and the omission is recorded in `strategy.yml` and the report.

With no GPU this falls back to sharded BRAKER and warns that the run may not
finish.

### Tools that do not apply here

*P. densiflora* is a **gymnosperm**. Two tools that look applicable are not:

- **Tiberius** ships Angiosperms (Mesangiospermae); no gymnosperm model.
- **NCBI EGAPx** supports Magnoliopsida; gymnosperms are out of scope.

Helixer's `land_plant` model is the only deep-learning track that generalises
here. This is recorded in `docs/clade_notes.md` and there is deliberately no
parameter for either tool.

## Practical notes

- Put `NXF_WORK` on the large filesystem. The launcher defaults to
  `/data/db/eukannot/work`; budget at least 3× genome size.
- `--publish_mode link` needs the work dir and the project output
  (`output/<project>/`) on one filesystem.
- Expect `GENOME_STATS` to take tens of minutes: it streams 21.7 Gb once.
  Composition (N%, softmask%) is measured on a 200 Mb prefix by default;
  length accounting is always exact.
- Infernal, Qualimap and TEsorter are auto-disabled above 2 Gb.
- Plan before you launch: `-entry strategy` prints the whole plan and runs
  nothing.
