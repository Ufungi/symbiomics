//
// HA<->HB assembly alignment and synteny/SNP calling -- the shared
// infrastructure behind both the diploid Salmon EM split (docs/haplotypes.md
// option 2) and the phASER allele-specific-expression route (option 3).
//
// minimap2 -x asm5/asm10/asm20 aligns two haplotype-resolved assemblies of the
// SAME individual; SyRI turns that alignment into synteny blocks and SNPs.
// Phase is known by construction here -- HA and HB are already the two
// separated haplotypes, so there is no statistical phasing step the way there
// would be calling heterozygous sites against one diploid reference.
//

process MINIMAP2_ASM {
    label 'process_high'
    tag "${ha.name} vs ${hb.name}"

    input:
    path ha
    path hb
    val  preset   // asm5 | asm10 | asm20 -- divergence-dependent, see docs/allele_specific_expression.md

    output:
    path 'ha_vs_hb.paf', emit: paf
    path 'versions.yml', emit: versions

    script:
    // Verified command shape (SyRI's own documented recipe): minimap2 without
    // -a/-ax already emits PAF by default; --eqx distinguishes matches from
    // mismatches in the CIGAR, which SyRI's SNP calling needs; SyRI then reads
    // this PAF directly with -F P.
    """
    minimap2 -x ${preset} --eqx -t ${task.cpus} '${ha}' '${hb}' > ha_vs_hb.paf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """

    stub:
    """
    touch ha_vs_hb.paf
    echo '"${task.process}": {minimap2: stub}' > versions.yml
    """
}

process SYRI {
    label 'process_high_memory'
    label 'process_long'
    tag "${ha.name} vs ${hb.name}"

    input:
    path paf
    path ha
    path hb

    output:
    path 'syri.out'     , emit: syri_out
    path 'syri.vcf'     , emit: vcf, optional: true
    path 'versions.yml' , emit: versions

    script:
    // -c/-r/-q/-F verified against SyRI's documented minimap2+PAF recipe.
    // --nc (threads) / --dir (output dir) are standard SyRI flags from general
    // documentation but not re-verified against this exact version this
    // session -- if either is renamed, this fails fast (unrecognised flag),
    // not silently.
    """
    syri -c '${paf}' -r '${ha}' -q '${hb}' -F P --nc ${task.cpus} --dir \$PWD

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        syri: \$(syri --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
    END_VERSIONS
    """

    stub:
    """
    printf 'chr1a\\t100\\t100\\tA\\tG\\tchr1b\\t100\\t100\\tSNP1\\t-\\tSNP\\t-\\n' > syri.out
    echo '"${task.process}": {syri: stub}' > versions.yml
    """
}
