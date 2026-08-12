//
// Allele-specific expression from the HA-only HISAT2 BAM this pipeline
// already produces -- no new alignment mode. See docs/allele_specific_expression.md
// for why standalone WASP (not STAR+WASP) and where the phased VCF comes from.
//
//   HA BAM (RNASEQ_ALIGN, unmodified)
//        |
//        v
//   WASP_FIND_INTERSECTING_SNPS  -> keep.bam (no SNP overlap, pass through)
//        |                          to.remap.bam + remap.fastq
//        v
//   HISAT2_ALIGN (reused, same index)   <- re-align the allele-flipped reads
//        |
//        v
//   WASP_FILTER_REMAPPED_READS  -> reads that remapped to the same place
//        |
//        v
//   WASP_MERGE_KEEP  -> bias-corrected BAM
//        |
//        v
//   PHASER (+ phased VCF)  -> per-variant haplotypic counts
//        |
//        v
//   PHASER_GENE_AE (+ HA gene BED)  -> per-gene allelic expression
//

include { VCF_BGZIP_INDEX; VCF_TO_WASP_SNPDIR                        } from '../../../modules/local/ase'
include { WASP_FIND_INTERSECTING_SNPS; WASP_FILTER_REMAPPED_READS    } from '../../../modules/local/ase'
include { WASP_MERGE_KEEP; PHASER; PHASER_GENE_AE; GFF3_TO_GENE_BED  } from '../../../modules/local/ase'
include { HISAT2_ALIGN as HISAT2_REMAP                               } from '../../../modules/local/align'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_REMAP                       } from '../../../modules/local/align'

workflow ALLELE_SPECIFIC_EXPRESSION {
    take:
    bam            // [ meta, bam, index ] -- the existing HA-only HISAT2 output
    phased_vcf     // path
    gff3_ha        // path
    hisat2_index   // path -- the SAME index the original alignment used
    index_prefix   // val
    use_mmap       // val

    main:
    ch_versions = Channel.empty()

    VCF_BGZIP_INDEX(phased_vcf)
    ch_versions = ch_versions.mix(VCF_BGZIP_INDEX.out.versions)

    VCF_TO_WASP_SNPDIR(phased_vcf)
    ch_versions = ch_versions.mix(VCF_TO_WASP_SNPDIR.out.versions)
    ch_snp_dir  = VCF_TO_WASP_SNPDIR.out.snp_dir.first()

    GFF3_TO_GENE_BED(gff3_ha)
    ch_versions = ch_versions.mix(GFF3_TO_GENE_BED.out.versions)

    WASP_FIND_INTERSECTING_SNPS(bam, ch_snp_dir)
    ch_versions = ch_versions.mix(WASP_FIND_INTERSECTING_SNPS.out.versions.first())

    // Re-align the allele-flipped reads with the SAME aligner and index the
    // original alignment used -- this IS the WASP bias check.
    HISAT2_REMAP(WASP_FIND_INTERSECTING_SNPS.out.remap_fastq, hisat2_index, index_prefix, use_mmap)
    ch_versions = ch_versions.mix(HISAT2_REMAP.out.versions.first())

    SAMTOOLS_SORT_REMAP(HISAT2_REMAP.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_REMAP.out.versions.first())

    to_remap_and_remapped = WASP_FIND_INTERSECTING_SNPS.out.to_remap_bam
        .join(SAMTOOLS_SORT_REMAP.out.bam, by: 0)

    WASP_FILTER_REMAPPED_READS(to_remap_and_remapped)
    ch_versions = ch_versions.mix(WASP_FILTER_REMAPPED_READS.out.versions.first())

    keep_and_filtered = WASP_FIND_INTERSECTING_SNPS.out.keep_bam
        .join(WASP_FILTER_REMAPPED_READS.out.bam, by: 0)

    WASP_MERGE_KEEP(keep_and_filtered)
    ch_versions = ch_versions.mix(WASP_MERGE_KEEP.out.versions.first())

    PHASER(WASP_MERGE_KEEP.out.bam, VCF_BGZIP_INDEX.out.vcf_gz.first(), VCF_BGZIP_INDEX.out.tbi.first())
    ch_versions = ch_versions.mix(PHASER.out.versions.first())

    PHASER_GENE_AE(PHASER.out.haplotypic_counts, GFF3_TO_GENE_BED.out.bed.first())
    ch_versions = ch_versions.mix(PHASER_GENE_AE.out.versions.first())

    emit:
    gene_ae  = PHASER_GENE_AE.out.gene_ae
    versions = ch_versions
}
