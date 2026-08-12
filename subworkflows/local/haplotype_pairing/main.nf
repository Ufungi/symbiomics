//
// HA<->HB pairing: one alignment step, two consumers.
//   - a phased VCF for phASER (docs/allele_specific_expression.md)
//   - a gene-level allele-pairing table for the diploid Salmon EM split
//     (docs/haplotypes.md option 2)
//
// Phase is known by construction (HA and HB are already the two separated
// haplotypes of one individual), so this is a deterministic assembly
// alignment, not a statistical phasing step.
//

include { MINIMAP2_ASM; SYRI } from '../../../modules/local/haplotype_align'

process SYRI_TO_PHASED_VCF {
    label 'process_single'
    tag "${syri_out.name}"

    input:
    path syri_out
    path genes_gff3   // NO_FILE to skip the gene-overlap restriction
    val  sample_name

    output:
    path 'ha_hb.phased.vcf', emit: vcf
    path 'versions.yml'    , emit: versions

    script:
    def genes = genes_gff3.name != 'NO_FILE' ? "--genes-gff3 '${genes_gff3}'" : ''
    """
    syri_to_phased_vcf.py --syri-out '${syri_out}' ${genes} \\
        --sample-name '${sample_name}' --out ha_hb.phased.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ha_hb.phased.vcf
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}

process SYRI_TO_ALLELE_TABLE {
    label 'process_single'
    tag "${syri_out.name}"

    input:
    path syri_out
    path ha_gff3
    path hb_gff3

    output:
    path 'ha_hb.allele_pairs.tsv', emit: table
    path 'versions.yml'          , emit: versions

    script:
    """
    syri_to_allele_table.py --syri-out '${syri_out}' \\
        --ha-gff3 '${ha_gff3}' --hb-gff3 '${hb_gff3}' \\
        --out ha_hb.allele_pairs.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ha_hb.allele_pairs.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}

workflow HAPLOTYPE_PAIRING {
    take:
    genome_ha    // path
    genome_hb    // path
    gff3_ha      // path
    gff3_hb      // path
    minimap_preset // val: asm5 | asm10 | asm20

    main:
    ch_versions = Channel.empty()

    MINIMAP2_ASM(genome_ha, genome_hb, minimap_preset)
    ch_versions = ch_versions.mix(MINIMAP2_ASM.out.versions)

    SYRI(MINIMAP2_ASM.out.paf, genome_ha, genome_hb)
    ch_versions = ch_versions.mix(SYRI.out.versions)

    SYRI_TO_PHASED_VCF(SYRI.out.syri_out, gff3_ha, 'HA_HB')
    ch_versions = ch_versions.mix(SYRI_TO_PHASED_VCF.out.versions)

    SYRI_TO_ALLELE_TABLE(SYRI.out.syri_out, gff3_ha, gff3_hb)
    ch_versions = ch_versions.mix(SYRI_TO_ALLELE_TABLE.out.versions)

    emit:
    syri_out     = SYRI.out.syri_out
    phased_vcf   = SYRI_TO_PHASED_VCF.out.vcf
    allele_table = SYRI_TO_ALLELE_TABLE.out.table
    versions     = ch_versions
}
