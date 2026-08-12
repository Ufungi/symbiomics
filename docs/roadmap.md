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
  dbCAN — all three against databases already on this lab's disk. **Verified**:
  a real run against real Swiss-Prot sequences and the real local eggNOG/CAZy
  databases produced correct, biologically sensible descriptions (including a
  genuine eggNOG ortholog hit InterProScan/Swiss-Prot alone would have missed).
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
