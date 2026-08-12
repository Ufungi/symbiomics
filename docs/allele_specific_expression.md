# Allele-specific expression

Two independent routes to the same question — is one haplotype's allele of a
gene expressed more than the other's — sharing one piece of infrastructure:
the HA↔HB assembly alignment (`docs/haplotypes.md` calls this "the pairing
step").

| route | mechanism | granularity |
|---|---|---|
| **A. Diploid Salmon EM, unsummed** | `docs/haplotypes.md` option 2, kept per-haplotype instead of summed | gene-level |
| **B. phASER** | WASP-filtered HISAT2 reads + read-backed haplotype counting | variant/gene-level |

Route A is already built (`--quant_engine salmon --transcript_fasta HA.CDS.fa,HB.CDS.fa`,
see `docs/haplotypes.md`). This page is about route B, and about the pairing
step both routes depend on.

## The pairing step (`-entry pairing`)

```bash
scripts/eukannot run . -entry pairing -profile singularity,local64 \
    --genome HA.genome.fasta --genome_hb HB.genome.fasta \
    --gff3_ha HA.gff3 --gff3_hb HB.gff3 \
    --minimap_preset asm10
```

```
HA.genome.fasta  ┐
                 ├─► minimap2 -x asm10 --eqx ─► ha_vs_hb.paf
HB.genome.fasta  ┘
                        │
                        ▼
                 SyRI -c paf -r HA -q HB -F P
                        │
          ┌─────────────┴──────────────┐
          ▼                             ▼
   syri.out SNP rows              syri.out SYN blocks
          │                             │
          ▼                             ▼
   phased VCF (GT=0|1,             gene-pairing table
   HA=REF, HB=ALT,                 (HA gene id <-> HB gene id,
   phase known by                  collinear within each block)
   construction)
```

**Phase is known by construction, not statistically inferred.** HA and HB are
already the two separated haplotypes of one individual — there is no
population-scale phasing uncertainty to resolve the way there would be if
heterozygous sites were called against one diploid reference. Every SNP from
the alignment gets `GT=0|1`, always, never `0/1`.

`--minimap_preset`: `asm5` (<1% divergence), `asm10` (a few %, the default —
appropriate for two haplotypes of one individual), `asm20` (up to ~10%,
cross-species territory).

**Gene pairing is positional within each SyRI syntenic (SYN) block**: the Nth
HA gene overlapping a block pairs with the Nth HB gene overlapping the same
block. This assumes genes appear in the same relative order on both sides of
a block, which is what "syntenic" means — reasonable for two haplotypes of one
individual, much less reliable across species. A block whose HA/HB gene counts
disagree leaves its genes unpaired rather than guessing which one is missing.

### Route B: only the phased VCF matters here

`--ase_vcf_source assembly` (default, the pairing step above) or
`--ase_vcf_source resequencing`.

**Why a second route exists**: the *P. densiflora* dataset ships its own
resequencing-derived genotypes
(`Variation_information_of_P.densiflora_accessions.txt.gz`). Where the RNA-seq
sample's accession is one the assembly authors resequenced, using their
published genotypes instead of this pipeline's own HA-vs-HB SNP calls is lower
review risk — that data is already peer-reviewed, not something this pipeline
computed and has to defend.

```bash
scripts/eukannot run . -entry variation \
    --resequencing_variation_file Variation_information_of_P.densiflora_accessions.txt.gz \
    --resequencing_accession <name>
```

**Its exact format was not known when this was written** — only the Figshare
description ("genotype information generated from resequencing analysis").
`bin/inspect_variation_file.py` follows the same discipline as v0.1's
`inspect_annotation.py`: read it, report the detected columns, and either
convert or fail with a specific message — never assume a layout that hasn't
been confirmed. If the file turns out to already be VCF, it passes straight
through. Route A (assembly-derived) does not depend on this file at all and
remains available regardless.

## The ASE workflow

```bash
scripts/eukannot run . -profile singularity,local64 \
    --input samplesheet.tsv --genome HA.genome.fasta --hisat2_index HA.genome \
    --run_ase true --ase_phased_vcf results_pairing/70_haplotype_pairing/ha_hb.phased.vcf \
    --gff3_ha HA.gff3 --outdir results
```

Runs downstream of the **existing, unmodified HA-only HISAT2 BAM** — this is
not a new alignment mode, just new consumers of the BAM v0.1 already produces.

```
HA BAM (RNASEQ_ALIGN, unmodified)
        │
        ▼
WASP_FIND_INTERSECTING_SNPS   -- flips the allele at each read's overlapping
        │                        heterozygous site, emits reads needing
        │                        re-alignment plus a pass-through set with no
        │                        overlapping SNP at all
        ▼
HISAT2 (reused, same index)   -- re-align the flipped reads
        ▼
WASP_FILTER_REMAPPED_READS    -- keep only reads that remap to the same place
        ▼
merge with the pass-through set, sort, index
        ▼
phASER (+ phased VCF)          -- per-variant haplotypic counts
        ▼
phASER Gene AE (+ HA gene BED) -- per-gene allelic expression
```

### Why standalone WASP, not STAR+WASP

STAR+WASP (`--waspOutputMode`) is the common recipe in the literature. This
pipeline's policy engine already refuses STAR above `--star_max_genome`
(2 Gb, `docs/decisions.md`) because its index is not buildable at conifer
scale — that constraint doesn't change for ASE. WASP's **original
implementation is aligner-agnostic**: `find_intersecting_snps.py` /
`filter_remapped_reads.py` work with any aligner, because the bias check is
"does this read map back to the same place after its allele is flipped and
it's re-aligned" — not anything STAR-specific. Re-aligning with HISAT2 (the
same index the original alignment used) is exactly equivalent to what
`--waspOutputMode` does internally, run as explicit steps.

## What is verified, and what is not, this session

Built and tested this session:
- `syri_to_phased_vcf.py` / `syri_to_allele_table.py`: unit-tested against a
  synthetic `syri.out` (correct SNP extraction with `0|1` GT, correct
  gene-overlap filtering, correct positional pairing within a block).
- `vcf_to_wasp_snp_dir.py`: unit-tested against a synthetic VCF.
- The full Nextflow DAG (`-entry pairing`, and `--run_ase true` wired into the
  main workflow): `-stub-run` verified end to end, including the real
  three-way BAM fork feeding both counting and the ASE branch simultaneously.
- Container tags for minimap2 (2.31), SyRI (1.8.2), and phASER (0.1.1ad5f89,
  bioconda) were individually confirmed to resolve.

**Not run against real data this session** (would need a multi-hour-to-day
whole-genome SyRI run, which was out of scope for this pass — the *P.
densiflora* HB assembly was downloaded but a full HA-vs-HB alignment was not
attempted):
- `minimap2`/`SyRI` on the actual 21.7 Gb HA/HB assemblies.
- The WASP mapping scripts (`find_intersecting_snps.py`,
  `filter_remapped_reads.py`) against real data. **They are not on
  bioconda or PyPI** (confirmed this session) — their dependencies (numpy,
  scipy, pysam) are pinned in `envs/wasp.yml`, but the scripts themselves
  still need to be staged into whatever environment runs them (`git clone
  https://github.com/bmvdgeijn/WASP.git`, add `mapping/` to `PATH`) or into a
  purpose-built container. **This is a real gap to close before route B runs
  for real** — see the risk table below.
- phASER's and `phaser_gene_ae`'s exact flag names beyond what was directly
  confirmed (`--bam`, `--vcf`, `--sample`, `--paired_end`, `--mapq`,
  `--baseq`, `--o` for `phaser.py`; `phaser_gene_ae.py`'s flags were not
  separately confirmed and are implemented from general documented usage).

## Risks

| # | risk | mitigation |
|---|---|---|
| 1 | WASP's mapping scripts have no package -- a container/env is pinned to their dependencies but the scripts themselves are not yet staged anywhere runnable | build a small image (`git clone` + the `pysam.yml`-equivalent deps) before running route B for real; tracked as an open gap, not silently assumed |
| 2 | phASER is old (2015-era), single bioconda build, py27 -- not verified for currency this session | pin the exact confirmed build (`0.1.1ad5f89--py27pl5321h9f5acd7_0`); if it breaks, GATK `ASEReadCounter` is the documented fallback for variant-level counting (loses phASER's read-backed haplotype phasing, but the phase is already known from the assembly alignment here, so the loss is smaller than in the typical population-ASE use case) |
| 3 | SyRI on two 21.7 Gb assemblies: runtime/memory not measured | budget it like every other huge-genome stage -- a time budget, and a documented fallback (route B via the resequencing file only, skipping the assembly-alignment route entirely) |
| 4 | Positional gene pairing within SYN blocks assumes collinearity | holds well for two haplotypes of one individual; a block with mismatched HA/HB gene counts is reported as unpaired rather than guessed, and a high unpaired fraction is flagged (`split_haplotype_counts.py`'s 30% warning threshold) |
| 5 | `phaser_gene_ae.py`'s exact flags unconfirmed | fails fast (unrecognised argument) rather than silently miscomputing if wrong; verify against `phaser_gene_ae.py --help` before a real run |
