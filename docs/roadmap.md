# Roadmap

Milestones are sized so something real runs early, and so the first release
already replaces the shell pipeline it supersedes.

## Naming: symbiomics, not eukannot

Renamed 2026-08-13 (`eukannot` → `symbiomics`, all files/paths/env vars —
`/data/db/eukannot`, the one pre-existing real database cache directory on
disk, was deliberately left unrenamed; see the git log for the mechanical
rename). The old name described what the pipeline *does* (eukaryote genome
annotation). The new name describes what it's *for*: this lab's actual
research subject is mycorrhizal symbiosis — *Pinus densiflora* (host) and
*Tricholoma matsutake* (ectomycorrhizal symbiont) — and "Symbiomics" is
meant as the umbrella name for that whole research direction, not just this
one pipeline.

**Scope honesty**: today, this pipeline still only does exactly what it did
under the old name — eukaryotic (fungal/plant) genome annotation and RNA
quantification, one organism's genome at a time. Simultaneous
fungi+bacteria community analysis and metabolomics — the rest of what
"Symbiomics" implies — have no design and no code yet. They belong to the
umbrella, not to this repo's current milestones. See "Beyond this pipeline"
at the bottom.

## v0.1.0 — skeleton + mRNA-seq arm ✅

- repo, launcher (`JAVA_CMD`), preflight, container/conda profiles
- samplesheet parsing, SE/PE auto-detection, resolved-sheet output
- **the policy engine** — `strategy.yml`, `-entry strategy`, budgets, fallbacks
- genome stats, size classes, clade exclusions
- fastp → strandedness inference → HISAT2 (`--large-index`, `--mm`, index reuse)
  → coordinate sort → CSI/BAI
- StringTie assemble + merge
- featureCounts / HTSeq counting, merged matrices, MultiQC

Two real bugs from the previous pipeline are fixed structurally rather than
patched: HTSeq is always given `-r pos` (its default `-r name` on a
coordinate-sorted BAM silently drops PE mates), and counting is over `exon`
features grouped by `gene_id`, never over `gene` spans.

## v0.2 — quantification engines, functional annotation, haplotype-aware analysis ✅

Reprioritized ahead of structural annotation: a published, haplotype-resolved
*Pinus densiflora* assembly (figshare 25546534, two full haplotypes, CC0)
already has good gene models and needs functional annotation and quantitative
analysis, not another round of gene prediction. See
`docs/functional_annotation.md`, `docs/haplotypes.md`, and
`docs/allele_specific_expression.md` for the full designs; this section is the
summary and what is/isn't verified.

**Salmon as a second quantification engine** (`--quant_engine salmon`) —
decoy-aware, bias-corrected, genome-alignment-free. HISAT2 stays required for
anything needing real alignments (phASER/WASP, StringTie, future BRAKER
evidence); Salmon is for fast bulk quantification and the diploid EM route.
Kallisto was compared and not added — it adds no capability Salmon lacks for
this pipeline's needs (no decoy-awareness, weaker bias modelling). **Verified**:
real end-to-end run (exact recovery of planted read depths); the disjoint-ID
guard; the GFF3-derived and identity-fallback gene-mapping ladder.

**`-entry functional`** — protein-in functional annotation, fully decoupled
from structural annotation. Two phases:
- **Phase A** (shipped, zero new downloads): DIAMOND-vs-Swiss-Prot, eggNOG-mapper,
  dbCAN v3 (HMMER + DIAMOND + eCAMI, consensus `#ofTools` column) — all three
  against databases already on this lab's disk. **Verified**: a real run
  against real Swiss-Prot sequences and the real local eggNOG/CAZy databases
  produced correct, biologically sensible descriptions (including a genuine
  eggNOG ortholog hit InterProScan/Swiss-Prot alone would have missed).
- **Phase B** (not yet wired): InterProScan 6, KofamScan, dbCAN v5, funannotate2
  + `funannotate2-addons` for the fungal branch (which already wraps eggNOG and
  InterProScan itself — not duplicated). Needs ~80 GB of new database
  provisioning.
- No AHRD anywhere: product naming is a direct best-hit ladder
  (`bin/assign_products.py`), not a weighted-blast synthesis step.

**Haplotype-resolved genome support** — `docs/haplotypes.md`'s three options,
now built:
1. Single-haplotype quantification (v0.1, unchanged).
2. Diploid Salmon EM (`--transcript_fasta HA.CDS.fa,HB.CDS.fa`) — summed
   (PAV-inclusive, reference-bias-free) and per-haplotype matrices.
   **Verified**: real run, correct arithmetic, and a real edge case (Salmon
   deduplicating byte-identical alleles at index time) found and handled
   correctly.
3. phASER allele-specific expression on WASP-filtered HISAT2 reads (not
   STAR+WASP — STAR is blocked above 2 Gb genomes in this pipeline's own
   policy). The phased VCF comes from either HA↔HB assembly alignment
   (minimap2 + SyRI, phase known by construction) or the dataset's own
   published resequencing genotypes (lower review risk, format not yet
   confirmed — `-entry variation` detects and reports rather than guessing).

**Why this stays justified rather than reinvented**: the closest published
end-to-end tool, ASET (2025, BMC Bioinformatics — a Nextflow ASE pipeline),
works from a single reference genome plus heterozygous SNP calls and
explicitly does not handle presence-absence variation genes. This pipeline's
core case — two physically assembled haplotypes, phase known by
genome-genome alignment rather than statistical phasing, PAV-aware summed
quantification — has no existing off-the-shelf equivalent. This is the part
of "Symbiomics" that is actually novel; the rest (aligning reads, calling
CAZymes, best-hit product naming) is deliberately unoriginal plumbing in
service of it.

**Verification honesty for this milestone**: the DAG (`-entry pairing`,
`-entry variation`, `--run_ase true`) is `-stub-run` verified end to end, and
every pure-Python script (`syri_to_phased_vcf.py`, `syri_to_allele_table.py`,
`vcf_to_wasp_snp_dir.py`, `split_haplotype_counts.py`,
`inspect_variation_file.py`) is unit-tested against synthetic data. **Not run
against real data this session**: minimap2/SyRI on the actual 21.7 Gb HA/HB
assemblies (a multi-hour-to-day computation), and the WASP mapping scripts
end-to-end (they are not packaged on bioconda/PyPI — confirmed this session —
and need a container built before route B is real-data-ready; tracked as an
open gap in `docs/allele_specific_expression.md`, not silently assumed done).

## v0.2.5 — `--organism` preset system

Not yet built. `--taxon` already exists but only branches one entry point
(`-entry functional`: `fungi` → funannotate2, `plant`/`animal`/`other` →
the shared core). The idea, raised and agreed 2026-08-13: extend that same
branching to *every* entry point through a single `--organism` flag that
picks sane per-organism defaults — which functional modules run, which
quant engine, and (once v0.3+ exists) which structural predictor track —
instead of `--taxon` covering only the one corner it covers today.

**Hard constraint carried over from this pipeline's existing decision-driven
design**: a preset sets *defaults*, never a lock. Every individual flag
(`--functional_modules`, `--quant_engine`, etc.) must remain independently
overridable after `--organism` picks its starting point — otherwise "this
species can only run this combination" silently replaces "this species
defaults to this combination," which breaks the policy-engine philosophy
already established in `strategy.yml`.

Not sized/scheduled in detail yet — this entry is a placeholder to track the
decision until it's actually designed.

## v0.3 — repeats + BRAKER (structural annotation)

Deferred behind v0.2 because it is not needed for the immediate *P.
densiflora* analysis (a good external annotation already exists) — but still
the next milestone for projects that need gene prediction from scratch.

- `repeats` subworkflow: `premasked`, `repeat_lib`, `earlgrey`, `edta`,
  `repeatmodeler_subsample` (+ `sample_contigs.py`, `salvage_repeatmodeler.py`),
  `red`, sharded RepeatMasker, masked-fraction guard, TEsorter
- `evidence_protein`: OrthoDB v12 partition resolution, `ncbi:` proteomes via
  gffread, miniprot
- BRAKER native, evidence-mode auto-derivation, TE/length/internal-stop filters
- `gff_finalize.py` — locus tags, sort, validate
- annotation QC: BUSCO, compleasm, OMArk, AGAT stats, gffcompare
- the regression harness against a known reference annotation

**Why BRAKER instead of MAKER**: verified 2026-08-13 — BRAKER3 outperforms
MAKER2 (and funannotate, and FINDER) by ~20 percentage points average
transcript-level F1 in its own published benchmark, with the gap *widest on
large, complex genomes* — directly relevant to this pipeline's 21.7 Gb
conifer target. See the "Deliberately deferred" table for the citation and
for a newer pipeline (Eukan) that was checked and not adopted.

Gate: reproduce a reference annotation with `--braker_version 3` first, then
measure and record the BRAKER 4 delta before promoting 4 to the default.

## v0.4 — multiple predictor tracks, consensus, functional Phase B

- Helixer (GPU), GALBA2, funannotate2 predict, GeMoMa
- TSEBRA / Mikado consensus, `consensus_method auto`
- **track budgets and graceful degradation wired to real processes**
- Functional Phase B: InterProScan 6, KofamScan, dbCAN v5, DeepSig, TMbed,
  Foldseek rescue tier, funannotate2 + `funannotate2-addons` fungal branch

## v0.5 — clade branches

- fungal: funannotate2 annotate, antiSMASH, MEROPS
- plant: Resistify (NLR), iTAK (TF/kinase), CYP450, OrthoFinder symbol transfer,
  Mercator4 import (a documented manual handoff, not a pretend pipeline step)

## v0.6 — huge genomes

- `make_shards.py`, sharded BRAKER train/predict/merge, sharded miniprot
- Helixer at 20+ Gb with scratch management
- end-to-end structural annotation on the 21.7 Gb conifer target

## v1.0 — polish

ncRNA, multi-genome batch via `--genomes`, `download_dbs` completeness,
container digest CI, nf-test coverage, full docs, Zenodo DOI.

## Deliberately deferred

| | why |
|---|---|
| sharded BRAKER | belongs to v0.6; a 22 Gb genome cannot ship early and pretending otherwise would be dishonest |
| Mikado | current build is a release candidate |
| EVidenceModeler | the author declared it unmaintained in 2024; legacy flag only |
| funannotate2-addons | very young (a dozen commits) |
| Kallisto | compared against Salmon (v0.2); adds no capability Salmon lacks here |
| DeepLoc, SignalP 6 | license-gated; DeepSig is the redistributable default |
| DeepTMHMM | sends sequences to a cloud service by default; TMbed is more accurate and permissively licensed |
| Tiberius, EGAPx | no gymnosperm models — see `docs/clade_notes.md` |
| WASP container build | scripts are unpackaged; dependencies pinned (`envs/wasp.yml`) but the scripts themselves are not yet staged into a runnable image — see `docs/allele_specific_expression.md` |
| MAKER | superseded for v0.3's purpose — BRAKER3 beats MAKER2 by ~20pp average transcript-F1 in its own benchmark, gap widest on large/complex genomes; [BRAKER3, *Genome Research* 2024](https://genome.cshlp.org/content/early/2024/05/28/gr.278090.123.full.pdf) |
| Eukan | checked 2026-08-13 — new (NAR Genomics & Bioinformatics, 2026-03), benchmarked against MAKER/BRAKER/GeMoMa on 17 genomes with statistically indistinguishable accuracy; wins on not-having-outlier-failures, not on raw power. Validated only up to ~400 Mb (rice) — no evidence at conifer scale. Watch, don't adopt yet |

## Beyond this pipeline — the Symbiomics umbrella

Not roadmap items — not designed, not scheduled, listed only so the name
doesn't imply more than the code does:

- **Bacterial community analysis** alongside the host and the fungal
  symbiont (16S/metagenomic arm for the mycorrhizosphere microbiome).
- **Metabolomics integration** (LC-MS/GC-MS), tying host-symbiont gene
  expression to the actual metabolic exchange symbiosis is about.

If either of these gets built, it's a separate pipeline/repo under the same
Symbiomics umbrella, not a feature bolted onto this one — this repo stays
scoped to genome annotation and RNA quantification (see
`docs/allele_specific_expression.md`'s and this file's existing scope notes;
also decided 2026-08-13: no DESeq2/edgeR-style DE/ASE significance calling
in this repo either, for the same reason).
