# mRNA-seq against a haplotype-resolved (diploid) assembly

Reference case: the published *Pinus densiflora* v1.0 assembly, which ships two
full haplotypes — **HA** and **HB** — each with its own genome FASTA, GFF3, CDS
and protein FASTA, ~44,233 genes apiece.

A diploid reference breaks the assumption every short-read counter is built on:
that a read belongs to one place. Almost every gene exists twice, at ~99%
identity, so almost every read maps twice.

## The trap

**Concatenating HA + HB into one FASTA and running HISAT2 → featureCounts/HTSeq
produces meaningless counts.** In detail:

- Every read from a conserved region aligns equally well to both haplotypes.
- HISAT2 reports it as multi-mapping and assigns a low MAPQ.
- featureCounts and HTSeq **discard multi-mapping reads by default**, so most of
  the library is thrown away.
- Turning multi-mapping on (`-M`) does the opposite and **double-counts** every
  allele.
- What survives is decided by arbitrary tie-breaking, not biology.

The failure is quiet: the pipeline finishes, the matrix looks normal, and the
assigned-read fraction is the only clue. Check it.

Chromosome names are a second hazard: if both haplotypes name a chromosome
`chr1a`, concatenation silently merges two different sequences. Verify the names
are disjoint before combining anything.

## The three real options

### 1. One haplotype only — the default for differential expression

Align to HA, count on `HA.gff3`. Nothing else changes.

```bash
scripts/eukannot run . -profile singularity,local64 \
    --genome P.densiflora_v1.0_HA.genome.fasta \
    --hisat2_index /path/P.densiflora_v1.0_HA.genome \
    --quant_gff P.densiflora_v1.0_HA.gff3 \
    --genome_size_class huge --taxon plant --clade gymnosperm \
    --input samplesheet.tsv --outdir results
```

**What you accept:** reads carrying an HB-specific allele align slightly worse,
and a few fail — reference bias. **Why it is usually fine:** the bias is a
property of the reference, so it is the *same in every sample*. Comparing
conditions, it largely cancels. What you cannot see are genes present only in
HB (PAV genes), which simply have no coordinates to be counted at.

**Use this when** the question is "which genes change between conditions".

### 2. Both haplotypes, EM quantification — unbiased totals, and ASE

Do not align to a concatenated genome. Quantify against a **transcriptome** of
both haplotypes with a tool that resolves multi-mapping probabilistically —
Salmon or RSEM — using the shipped `HA.CDS.fa` and `HB.CDS.fa`.

```
cat HA.CDS.fa HB.CDS.fa > diploid.CDS.fa      # check IDs are disjoint first
salmon index -t diploid.CDS.fa -i diploid_idx
salmon quant -i diploid_idx -l A -1 R1.fq.gz -2 R2.fq.gz -o quant
```

Then either:
- **sum each allele pair** → total gene expression with no reference bias, and
  PAV genes included; or
- **keep the pair separate** → allele-specific expression, which is the whole
  point of a haplotype-resolved assembly.

This needs an **HA ↔ HB allele pairing table**, which the dataset does not ship.
Two ways to get one:

1. Check whether the gene IDs are already parallel. HA genes are named
   `Pd01G00010A`; if HB uses `Pd01G00010B` for the allele, pairing is free:
   ```bash
   grep -c "" <(awk -F'\t' '$3=="gene"' HB.gff3)      # same 44,233?
   # then pair on the ID stem, and confirm on a sample by position
   ```
   Do not assume it — confirm the pairing on coordinates for a few hundred genes
   before trusting it genome-wide.
2. Otherwise derive it: reciprocal best DIAMOND hits between `HA.PEP.fa` and
   `HB.PEP.fa`, restricted to syntenic positions.

**Use this when** you care about allelic imbalance, or when PAV genes matter.

### 3. Allele-specific expression from one reference

Align to HA, then count reads over heterozygous sites rather than over genes.
The dataset's `Variation_information_of_P.densiflora_accessions.txt.gz` gives
the variant set. Tools: GATK `ASEReadCounter` or phASER, with **WASP** to remove
the mapping bias that otherwise inflates the reference allele.

**Use this when** you want per-site allelic ratios in a specific accession.

## Recommendation

| goal | approach |
|---|---|
| DE between conditions | **HA only** (option 1) |
| total expression without reference bias, or PAV genes | Salmon on HA+HB CDS (option 2) |
| allelic imbalance | option 2 kept unsummed, or option 3 |
| anything | **never** HA+HB concatenated genome + featureCounts/HTSeq |

Start with option 1. It answers the usual question, it is immediately runnable,
and it is what the field does. Move to option 2 only when the biology you are
after is specifically allelic.

## Status in this pipeline

v0.1.0 implements option 1 directly — it is just a normal run against one
haplotype. Options 2 and 3 need a Salmon/RSEM quantification path and an allele
pairing step, which are not built yet; see `docs/roadmap.md`.

Until a guard exists, **check the assigned-read fraction** in
`95_quantify/counts_*.stats.tsv` and the alignment summaries in `40_align/log/`.
On a single haplotype, expect the usual rates. A sharp drop with a large
multi-mapping fraction means a duplicated reference slipped in.
