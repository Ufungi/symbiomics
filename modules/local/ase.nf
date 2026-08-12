//
// Allele-specific expression: WASP mapping-bias filtering (standalone, reusing
// HISAT2 for the remap step -- see docs/allele_specific_expression.md for why
// not STAR+WASP) followed by phASER.
//
// Verify-at-deployment note: WASP's exact per-run output filenames were
// confirmed from its documented example (PREFIX.keep.bam / PREFIX.to.remap.bam
// / PREFIX.remap.fq[12].gz inside --output_dir) but the prefixing convention
// itself was not exercised against a real run this session -- processes below
// glob for the documented suffixes rather than hardcoding an assumed prefix,
// so a different prefix convention still resolves correctly.
//

process VCF_BGZIP_INDEX {
    label 'process_single'
    tag "${vcf.name}"

    input:
    path vcf

    output:
    path '*.vcf.gz'    , emit: vcf_gz
    path '*.vcf.gz.tbi', emit: tbi
    path 'versions.yml', emit: versions

    script:
    """
    bgzip -c '${vcf}' > ${vcf.baseName}.vcf.gz
    tabix -p vcf ${vcf.baseName}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htslib: \$(bgzip --version | head -1 | grep -oE '[0-9]+\\.[0-9]+' | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch ${vcf.baseName}.vcf.gz ${vcf.baseName}.vcf.gz.tbi
    echo '"${task.process}": {htslib: stub}' > versions.yml
    """
}

process VCF_TO_WASP_SNPDIR {
    label 'process_single'
    tag "${vcf.name}"

    input:
    path vcf

    output:
    path 'wasp_snps', emit: snp_dir
    path 'versions.yml', emit: versions

    script:
    """
    vcf_to_wasp_snp_dir.py --vcf '${vcf}' --out-dir wasp_snps

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p wasp_snps && touch wasp_snps/chr1.snps.txt
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}

process WASP_FIND_INTERSECTING_SNPS {
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(index)
    path  snp_dir

    output:
    tuple val(meta), path('wasp_out/*.keep.bam')      , emit: keep_bam
    tuple val(meta), path('wasp_out/*.to.remap.bam')  , emit: to_remap_bam
    tuple val(meta), path('wasp_out/*.remap.fq*.gz')  , emit: remap_fastq
    path 'versions.yml', emit: versions

    script:
    def paired = meta.single_end ? '' : '--is_paired_end'
    """
    mkdir -p wasp_out
    find_intersecting_snps.py ${paired} --is_sorted \\
        --output_dir wasp_out --snp_dir '${snp_dir}' '${bam}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp: 2015-era mapping scripts (bmvdgeijn/WASP), not versioned by bioconda
    END_VERSIONS
    """

    stub:
    """
    mkdir -p wasp_out
    touch wasp_out/${meta.id}.keep.bam wasp_out/${meta.id}.to.remap.bam
    ${meta.single_end
        ? "touch wasp_out/${meta.id}.remap.fq.gz"
        : "touch wasp_out/${meta.id}.remap.fq1.gz wasp_out/${meta.id}.remap.fq2.gz"}
    echo '"${task.process}": {wasp: stub}' > versions.yml
    """
}

process WASP_FILTER_REMAPPED_READS {
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(to_remap_bam), path(remapped_bam)

    output:
    tuple val(meta), path('*.wasp_filtered.bam'), emit: bam
    path 'versions.yml', emit: versions

    script:
    """
    filter_remapped_reads.py '${to_remap_bam}' '${remapped_bam}' ${meta.id}.wasp_filtered.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wasp: 2015-era mapping scripts (bmvdgeijn/WASP), not versioned by bioconda
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.wasp_filtered.bam
    echo '"${task.process}": {wasp: stub}' > versions.yml
    """
}

process WASP_MERGE_KEEP {
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(keep_bam), path(filtered_bam)

    output:
    tuple val(meta), path('*.wasp.sorted.bam'), emit: bam
    path 'versions.yml', emit: versions

    script:
    """
    samtools merge -f -@ ${task.cpus} merged.bam '${keep_bam}' '${filtered_bam}'
    samtools sort -@ ${task.cpus} -o ${meta.id}.wasp.sorted.bam merged.bam
    samtools index -c ${meta.id}.wasp.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.wasp.sorted.bam
    echo '"${task.process}": {samtools: stub}' > versions.yml
    """
}

process PHASER {
    label 'process_medium'
    tag "${meta.id}"

    // Indexes the bam itself (below) rather than taking one as an input --
    // WASP_MERGE_KEEP already writes a CSI alongside it, but re-indexing here
    // keeps this process self-contained regardless of which index format its
    // upstream produced.
    input:
    tuple val(meta), path(bam)
    path  vcf_gz
    path  vcf_tbi

    output:
    tuple val(meta), path('*.haplotypic_counts.txt'), emit: haplotypic_counts
    path 'versions.yml', emit: versions

    script:
    def paired = meta.single_end ? 0 : 1
    """
    samtools index '${bam}'
    phaser.py --bam '${bam}' --vcf '${vcf_gz}' --sample '${params.ase_sample_name}' \\
        --paired_end ${paired} --mapq 20 --baseq 10 \\
        --threads ${task.cpus} --o ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phaser: 0.1.1ad5f89 (bioconda phaser, py27 -- see docs/allele_specific_expression.md)
    END_VERSIONS
    """

    stub:
    """
    printf 'contig\\tstart\\tstop\\tvariants\\thap_A_count\\thap_B_count\\n' > ${meta.id}.haplotypic_counts.txt
    echo '"${task.process}": {phaser: stub}' > versions.yml
    """
}

process PHASER_GENE_AE {
    label 'process_single'
    tag "${meta.id}"

    input:
    tuple val(meta), path(haplotypic_counts)
    path  gene_bed

    output:
    tuple val(meta), path('*.gene_ae.txt'), emit: gene_ae
    path 'versions.yml', emit: versions

    script:
    """
    phaser_gene_ae.py --haplotypic_counts '${haplotypic_counts}' \\
        --features '${gene_bed}' --o ${meta.id}.gene_ae.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phaser: 0.1.1ad5f89 (bioconda phaser)
    END_VERSIONS
    """

    stub:
    """
    printf 'contig\\tstart\\tstop\\tname\\taCount\\tbCount\\ttotalCount\\n' > ${meta.id}.gene_ae.txt
    echo '"${task.process}": {phaser: stub}' > versions.yml
    """
}

process GFF3_TO_GENE_BED {
    label 'process_single'
    tag "${gff3.name}"

    input:
    path gff3

    output:
    path 'genes.bed', emit: bed
    path 'versions.yml', emit: versions

    script:
    """
    awk -F'\\t' '\$0 !~ /^#/ && \$3 == "gene" {
        id = \$9; sub(/^.*ID=/, "", id); sub(/;.*\$/, "", id);
        print \$1"\\t"(\$4-1)"\\t"\$5"\\t"id
    }' '${gff3}' > genes.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: builtin
    END_VERSIONS
    """

    stub:
    """
    touch genes.bed
    echo '"${task.process}": {awk: stub}' > versions.yml
    """
}
