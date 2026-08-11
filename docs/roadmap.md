# Roadmap

Milestones are sized so something real runs early, and so the first release
already replaces the shell pipeline it supersedes.

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

## v0.2 — repeats + BRAKER

- `repeats` subworkflow: `premasked`, `repeat_lib`, `earlgrey`, `edta`,
  `repeatmodeler_subsample` (+ `sample_contigs.py`, `salvage_repeatmodeler.py`),
  `red`, sharded RepeatMasker, masked-fraction guard, TEsorter
- `evidence_protein`: OrthoDB v12 partition resolution, `ncbi:` proteomes via
  gffread, miniprot
- BRAKER native, evidence-mode auto-derivation, TE/length/internal-stop filters
- `gff_finalize.py` — locus tags, sort, validate
- annotation QC: BUSCO, compleasm, OMArk, AGAT stats, gffcompare
- the regression harness against a known reference annotation

Gate: reproduce a reference annotation with `--braker_version 3` first, then
measure and record the BRAKER 4 delta before promoting 4 to the default.

## v0.3 — multiple tracks, consensus, functional core

- Helixer (GPU), GALBA2, funannotate2 predict, GeMoMa
- TSEBRA / Mikado consensus, `consensus_method auto`
- **track budgets and graceful degradation wired to real processes**
- functional core: InterProScan 6, eggNOG-mapper, UPIMAPI (offline), DIAMOND vs
  Swiss-Prot, dbCAN, KofamScan, DeepSig, TMbed, **AHRD**, Foldseek rescue tier
- `assign_products.py` / `qc_products.py` — the ladder that keeps the final GFF3
  from being half "hypothetical protein"

## v0.4 — clade branches

- fungal: funannotate2 annotate, antiSMASH, MEROPS
- plant: Resistify (NLR), iTAK (TF/kinase), CYP450, OrthoFinder symbol transfer,
  Mercator4 import (a documented manual handoff, not a pretend pipeline step)

## v0.5 — huge genomes

- `make_shards.py`, sharded BRAKER train/predict/merge, sharded miniprot
- Helixer at 20+ Gb with scratch management
- end-to-end on the 21.7 Gb conifer target

## v1.0 — polish

ncRNA, multi-genome batch via `--genomes`, `download_dbs` completeness,
container digest CI, nf-test coverage, full docs, Zenodo DOI.

## Deliberately deferred

| | why |
|---|---|
| sharded BRAKER | belongs to v0.5; a 22 Gb genome cannot ship in v0.1 and pretending otherwise would be dishonest |
| Mikado | current build is a release candidate |
| EVidenceModeler | the author declared it unmaintained in 2024; legacy flag only |
| funannotate2-addons | very young (a dozen commits) |
| DeepLoc, SignalP 6 | license-gated; DeepSig is the redistributable default |
| DeepTMHMM | sends sequences to a cloud service by default; TMbed is more accurate and permissively licensed |
| Tiberius, EGAPx | no gymnosperm models — see `docs/clade_notes.md` |
