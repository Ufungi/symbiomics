<p align="center">
  <img src="https://img.shields.io/badge/platform-Linux-blue?logo=linux&logoColor=white" alt="Platform">
  <img src="https://img.shields.io/badge/nextflow-%E2%89%A524.10-23aa62?logo=nextflow&logoColor=white" alt="Nextflow">
  <img src="https://img.shields.io/badge/run%20with-singularity%20%7C%20docker%20%7C%20conda-blue?logo=singularity&logoColor=white" alt="Containers">
  <img src="https://img.shields.io/badge/language-Nextflow%20%7C%20Python%20%7C%20Bash-informational?logo=gnu-bash&logoColor=white" alt="Language">
  <img src="https://img.shields.io/badge/license-MIT-lightgrey" alt="License">
</p>

<p align="right"><strong>English</strong> | <a href="README.md">한국어</a></p>

# eukannot

A decision-driven Nextflow pipeline for eukaryote mRNA-seq quantification and
genome annotation. It reads the genome first, decides how to handle it, and
tells you why — before spending a single core-hour.

```bash
scripts/eukannot run . -profile singularity,local64 \
    --input samplesheet.tsv --genome genome.fasta --outdir results
```

Built and verified against two very different targets: *Tricholoma matsutake*
(161 Mb, softmasked fungal genome) and *Pinus densiflora* (21.74 Gb,
haplotype-resolved conifer genome, 2,006 contigs, longest contig 2.008 Gb). The
same command runs on both; the pipeline picks different tools underneath.

---

## ✨ Key Features

*   **A policy engine, not a tool list** — after measuring the genome,
    [`decide_strategy.py`](bin/decide_strategy.py) writes `strategy.yml`
    recording, for every choice, what will run, why, what was rejected, and
    which flag overrides it. `-entry strategy` resolves the whole plan and
    stops, so a multi-week run can be reviewed before it starts.
*   **Graceful degradation, not brittle failure** — heavy tracks carry wall-clock
    budgets. One that overruns is dropped, the run finishes with what is left,
    and the omission is written into the report instead of being silently lost.
*   **Align once** — one HISAT2 BAM channel forks to annotation evidence,
    StringTie assembly and read counting, so a 63-sample project never
    re-aligns the same reads twice.
*   **Two quantification engines** — [HISAT2](https://github.com/DaehwanKimLab/hisat2)
    for real alignments (needed for allele-specific expression and future
    structural-annotation evidence), or [Salmon](https://github.com/COMBINE-lab/salmon)
    for fast, bias-corrected, decoy-aware transcript quantification. Pick per run
    with `--quant_engine`.
*   **Structural annotation is optional, not assumed** — `-entry functional`
    takes a protein FASTA straight to functional annotation. If a good gene
    model already exists (as it does for *P. densiflora*), there is no need to
    run BRAKER just to get back to where you started.
*   **Annotation-free strandedness inference** — reads splice-junction motifs
    (`GT..AG` / `CT..AC`) directly out of the genome, because at this point in
    the pipeline there is no transcriptome yet to run Salmon- or RSeQC-style
    inference against.
*   **Haplotype-aware by design, not by accident** — a haplotype-resolved
    diploid assembly breaks the one-read-one-locus assumption every counter is
    built on. `docs/haplotypes.md` names the trap; the pipeline implements the
    three defensible ways around it, including a real allele-specific-expression
    workflow ([phASER](https://github.com/secastel/phaser) on
    [WASP](https://github.com/bmvdgeijn/WASP)-filtered reads).

---

## Pipeline overview

```
samplesheet.tsv ──► INPUT_CHECK (SE/PE auto-detect) ─────────────────┐
                                                                       │
genome.fasta ──► GENOME_PREP ──► DECIDE_STRATEGY ──► strategy.yml ────┤
                                                                       │
                        ┌──────────────────────────────────────────────┘
                        ▼
              fastp ─► strandedness probe ─► HISAT2 ─► sort ─► CSI/BAI index
                                                                       │
                 ┌─────────────────────────────┬─────────────────────┼──────────────┐
                 ▼                             ▼                     ▼              ▼
          StringTie assemble          featureCounts / HTSeq    (evidence,     phASER ASE
                                       or Salmon quant           v0.3+)        (optional,
                                                                                haplotype-
                                                                                resolved
                                                                                genomes)
                                              │
                                              ▼
                                        MultiQC report

  ── independent path, no genome required ──
  protein FASTA ──► -entry functional ──► taxon=fungi:  funannotate2 + f2a (not yet wired)
                                           taxon=other:  DIAMOND/Swiss-Prot + eggNOG-mapper +
                                                          dbCAN v3 (HMMER+DIAMOND+eCAMI)  [Phase A]
                                                          + InterProScan + KofamScan       [Phase B]
                                        ──► product-name ladder (no AHRD) ──► functional TSV
```

---

## Prerequisites

| | requirement |
|---|---|
| OS | Linux |
| Orchestrator | [Nextflow](https://www.nextflow.io/) ≥ 24.10 (needs Java 17–24) |
| Runtime | [Singularity](https://sylabs.io/singularity/) or Docker (recommended), or Conda (partial coverage — see `conf/conda.config`) |
| Storage | genome-dependent; budget ≥ 3× genome size for the work directory |
| GPU | optional — only used by later-milestone tracks (Helixer, TMbed) |

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/Ufungi/eukannot.git
cd eukannot
```

### 2. Check the environment

```bash
scripts/eukannot preflight
```

Verifies Java version, Nextflow, Singularity/Docker, GPU, disk space, and
database provisioning — each failure prints the exact fix. `scripts/eukannot`
is the supported entry point rather than calling `nextflow` directly: it
points `JAVA_CMD` at the JDK bundled with the Nextflow conda environment,
because Nextflow needs Java 17–24 and many hosts default to Java 11.

```bash
which nextflow    # should resolve once `conda activate nextflow` — or just
                   # use scripts/eukannot, which does this for you
```

### 3. (Optional) provision databases for functional annotation

```bash
scripts/eukannot run . -entry download_dbs --db_dir /data/db/eukannot
```

Not required for the mRNA-seq arm. Only needed for `-entry functional`'s
Phase B modules (InterProScan 6, dbCAN v5, KofamScan) — Phase A
(DIAMOND-vs-Swiss-Prot + eggNOG-mapper + dbCAN v3) runs against databases this
lab already has on disk. See `docs/functional_annotation.md`.

---

## Input files

| File | Required for | Description |
|---|---|---|
| `samplesheet.tsv` | mRNA-seq arm | One row per sample: `sample_id`, `fastq_1`/`fastq_2` or `bam`, optional `layout`/`strandedness`/`use_for`. SE/PE is auto-detected. See `assets/samplesheet.example.tsv`. |
| `genome.fasta` | mRNA-seq arm, structural annotation | Reference genome. Soft-masked or not — the pipeline masks it if needed and refuses to trust a `--premasked` claim with no lowercase in it. |
| `proteome.fasta` | `-entry functional` | Predicted protein sequences. No genome or samplesheet required for this path. |
| `genomes.tsv` | multi-genome batches | One row per project: fasta, taxon, masker, per-project overrides. See `assets/genomes.example.tsv`. |
| `reference_proteomes.tsv` | product-name ladder | Weighted, labelled reference proteomes (e.g. Swiss-Prot, a close relative) consumed by the functional-annotation product ladder. See `assets/reference_proteomes.example.tsv`. |

---

## Usage

**Full run — align, assemble, count:**

```bash
scripts/eukannot run . -profile singularity,local64 \
    --input samplesheet.tsv --genome genome.fasta \
    --taxon plant --outdir results
```

**Quantification only, with Salmon instead of HISAT2** (no BAM produced,
decoy-aware bias correction):

```bash
scripts/eukannot run . -profile singularity,local64 \
    --input samplesheet.tsv --quant_engine salmon \
    --transcript_fasta HA.CDS.fa --outdir results
```

**Functional annotation only — no structural annotation, no genome:**

```bash
scripts/eukannot run . -entry functional -profile singularity,local64 \
    --proteome HA.PEP.fa --genome_id Pinde_HA --taxon plant \
    --outdir results_functional
```

**Review the execution plan before committing to a multi-day run:**

```bash
scripts/eukannot run . -entry strategy \
    --genome genome.fasta --taxon plant --clade gymnosperm
```

More: `docs/usage.md` (parameters, output layout), `docs/decisions.md` (how
the policy engine chooses), `docs/haplotypes.md` (diploid/haplotype-resolved
genomes), `docs/functional_annotation.md`.

---

## Step reference

| Step | What it does | Runs when |
|---|---|---|
| `INPUT_CHECK` | samplesheet parsing, SE/PE inference, validation | always |
| `GENOME_PREP` / `DECIDE_STRATEGY` | genome stats, size class, `strategy.yml` | always (any entry that reads a genome) |
| `RNASEQ_ALIGN` | fastp, strandedness inference, HISAT2, sort, index | `--steps` includes `align` |
| `RNASEQ_ASSEMBLE` | StringTie per-sample assembly + merge | `assemble` |
| `QUANTIFY` | featureCounts / HTSeq / Salmon, merged matrices | `quantify` |
| `FUNCTIONAL` | DIAMOND-vs-Swiss-Prot / eggNOG-mapper / dbCAN v3 (HMMER+DIAMOND+eCAMI) — Phase A, shipped; InterProScan / KofamScan / funannotate2 — Phase B, not yet wired | `-entry functional` |
| `HAPLOTYPE_PAIRING` | minimap2 + SyRI: HA↔HB phased VCF and allele-pairing table | `-entry pairing`, standalone (its output feeds the two rows below) |
| `QUANTIFY` (Salmon, diploid mode) | summed + per-haplotype allele matrices | `--quant_engine salmon --transcript_fasta HA.CDS.fa,HB.CDS.fa --haplotype_pairs <pairing output>` |
| `ALLELE_SPECIFIC_EXPRESSION` | WASP-filtered HISAT2 reads → phASER Gene AE | `--run_ase true --ase_phased_vcf <pairing output>` |
| `MULTIQC` | consolidated report, `versions.yml`, `strategy.yml` | always |

---

## Tools used

| Tool | Role | Reference |
|---|---|---|
| [Nextflow](https://www.nextflow.io/) | workflow engine | Di Tommaso et al., *Nat. Biotechnol.* 2017 |
| [HISAT2](https://github.com/DaehwanKimLab/hisat2) | spliced short-read alignment | Kim et al., *Nat. Biotechnol.* 2019 |
| [Salmon](https://github.com/COMBINE-lab/salmon) | decoy-aware, bias-corrected transcript quantification | Patro et al., *Nat. Methods* 2017 |
| [StringTie](https://github.com/gpertea/stringtie) | transcript assembly | Pertea et al., *Nat. Biotechnol.* 2015 |
| [Subread/featureCounts](https://subread.sourceforge.net/) | exon-level read counting | Liao et al., *Bioinformatics* 2014 |
| [HTSeq](https://htseq.readthedocs.io/) | exon-level read counting (cross-check) | Anders et al., *Bioinformatics* 2015 |
| [samtools](https://www.htslib.org/) | BAM sort/index/stats | Danecek et al., *GigaScience* 2021 |
| [fastp](https://github.com/OpenGene/fastp) | read trimming/QC | Chen et al., *Bioinformatics* 2018 |
| [DIAMOND](https://github.com/bbuchfink/diamond) | Swiss-Prot best-hit product transfer (Phase A; UPIMAPI's own ID-mapping step is the Phase B upgrade to this slot) | Buchfink et al., *Nat. Methods* 2021 |
| [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) | orthology-based GO/KEGG/description transfer | Cantalapiedra et al., *Mol. Biol. Evol.* 2021 |
| [InterProScan](https://github.com/ebi-pf-team/interproscan6) | domain/family/GO annotation (Phase B) | Jones et al., *Bioinformatics* 2014 |
| [dbCAN / run_dbcan](https://github.com/bcb-unl/run_dbcan) (v3, HMMER+DIAMOND+eCAMI consensus) | CAZyme annotation | Zheng et al., *Nucleic Acids Res.* 2023 |
| [funannotate2](https://github.com/nextgenusfs/funannotate2) | fungal structural + functional annotation | Palmer & Stajich |
| [minimap2](https://github.com/lh3/minimap2) | assembly-to-assembly alignment (HA↔HB) | Li, *Bioinformatics* 2018 |
| [SyRI](https://github.com/schneebergerlab/syri) | synteny and structural-variant calling between two assemblies | Goel et al., *Genome Biol.* 2019 |
| [WASP](https://github.com/bmvdgeijn/WASP) | mapping-bias filtering for allele-specific reads | van de Geijn et al., *Nat. Methods* 2015 |
| [phASER](https://github.com/secastel/phaser) | read-backed haplotype phasing and allelic expression | Castel et al., *Nat. Commun.* 2016 |
| [MultiQC](https://multiqc.info/) | consolidated QC report | Ewels et al., *Bioinformatics* 2016 |

Full citation list with DOIs: `docs/citations.md` (filled in as each tool's
module is verified — see `scripts/check_containers.sh`).

---

## Status

**v0.1.0** — the mRNA-seq arm (samplesheet, strandedness inference, HISAT2,
StringTie, counting, the policy engine) is built and tagged.

**v0.2** — Salmon as a second quantification engine, `-entry functional`
(protein-in functional annotation, decoupled from structural annotation), and
haplotype-resolved-genome support: HA↔HB pairing (`-entry pairing`), diploid
Salmon quantification, and phASER allele-specific expression (`--run_ase`).
Salmon and functional-annotation Phase A (DIAMOND-vs-Swiss-Prot transfer,
eggNOG-mapper, dbCAN v3 with HMMER+DIAMOND+eCAMI) are verified against real
data this pass;
haplotype pairing and the ASE workflow are `-stub-run` verified and
unit-tested on synthetic data, with a known open gap (a WASP container still
needs building — see `docs/allele_specific_expression.md`) before route B runs
for real. See `docs/roadmap.md`.

Structural annotation (repeats, BRAKER, Helixer, consensus) is deliberately
later — a good external annotation already exists for this pipeline's
reference target, so functional annotation and quantification were
prioritized first.

---

## License

MIT. See `LICENSE`.

---

*eukaryote genome annotation · mRNA-seq · Nextflow · haplotype-resolved genomes · allele-specific expression · Pinus densiflora · Tricholoma matsutake*
