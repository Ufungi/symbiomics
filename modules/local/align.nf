//
// Alignment: HISAT2 -> coordinate-sorted CRAM -> index -> stats.
//
// Two things here are load-bearing and were bugs in the pipeline this replaces:
//   * SAM is never written to disk; hisat2 pipes straight into samtools sort.
//   * Alignments are stored as CRAM (reference-compressed), so the index is
//     always CRAI -- CRAM has no BAI 512 Mb ceiling.
//

process HISAT2_BUILD {
    label 'process_high'
    tag "${fasta.name}"

    input:
    path fasta
    val  large_index

    output:
    path 'hisat2_index'      , emit: index
    path 'versions.yml'      , emit: versions

    script:
    def large = large_index ? '--large-index' : ''
    """
    mkdir -p hisat2_index
    hisat2-build -p ${task.cpus} ${large} '${fasta}' hisat2_index/genome

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$(hisat2 --version | head -1 | sed 's/.*version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p hisat2_index && touch hisat2_index/genome.1.ht2
    echo '"${task.process}": {hisat2: stub}' > versions.yml
    """
}

process HISAT2_ALIGN {
    label 'process_high'
    tag "${meta.id}"

    // NOTE: Nextflow <=24 accepted a per-process containerOptions/runOptions
    // here to bind each task's work dir to /tmp (hisat2 FIFO collisions);
    // Nextflow 26 removed that directive, so this process runs without it.

    input:
    tuple val(meta), path(reads)
    path  index
    val   index_prefix
    val   use_mmap
    path  genome_fasta

    output:
    tuple val(meta), path('*.sorted.cram'), emit: bam
    tuple val(meta), path('*.summary.txt') , emit: summary
    path 'versions.yml'                    , emit: versions

    script:
    // --mm memory-maps the index so N concurrent aligners share one page-cache
    // copy. Without it, a 38 GB index times 13 forks exceeds any sane host.
    def mmap = use_mmap ? '--mm' : ''
    def extra = params.hisat2_extra ?: ''
    def strand = ''
    if( meta.strandedness == 'forward' ) strand = meta.single_end ? '--rna-strandness F' : '--rna-strandness FR'
    if( meta.strandedness == 'reverse' ) strand = meta.single_end ? '--rna-strandness R' : '--rna-strandness RF'

    def input_args = meta.single_end ? "-U ${reads.join(',')}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    hisat2 -p ${task.cpus} -x ${index}/${index_prefix} ${input_args} \\
        ${strand} ${mmap} ${extra} \\
        --new-summary --summary-file ${meta.id}.summary.txt \\
        | samtools sort -@ 2 -m ${params.sort_mem} -O cram --reference '${genome_fasta}' -o ${meta.id}.sorted.cram -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$(hisat2 --version | head -1 | sed 's/.*version //')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sorted.cram ${meta.id}.summary.txt
    echo '"${task.process}": {hisat2: stub}' > versions.yml
    """
}

process SAMTOOLS_SORT {
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam)
    path  genome_fasta

    output:
    tuple val(meta), path('*.sorted.cram'), emit: bam
    path 'versions.yml'                  , emit: versions

    script:
    """
    samtools sort -@ ${task.cpus} -m ${params.sort_mem} \\
        -O cram --reference '${genome_fasta}' \\
        -T ${meta.id}.tmp -o ${meta.id}.sorted.cram '${bam}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sorted.cram
    echo '"${task.process}": {samtools: stub}' > versions.yml
    """
}

process SAMTOOLS_INDEX {
    label 'process_low'
    tag "${meta.id}"

    input:
    tuple val(meta), path(cram)

    output:
    tuple val(meta), path(cram), path('*.crai'), emit: bam
    path 'versions.yml'                       , emit: versions

    script:
    // CRAM is always indexed with CRAI: it has no BAI 512 Mb ceiling, so the
    // csi/bai/both strategy decision is moot for CRAM.
    """
    samtools index -@ ${task.cpus} '${cram}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${cram}.crai
    echo '"${task.process}": {samtools: stub}' > versions.yml
    """
}

process BAM_STATS {
    label 'process_low'
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(index)
    path  genome_fasta

    output:
    path '*.stats'      , emit: stats
    path '*.flagstat'   , emit: flagstat
    path '*.idxstats'   , emit: idxstats
    path 'versions.yml' , emit: versions

    script:
    """
    samtools stats    --reference '${genome_fasta}' -@ ${task.cpus} '${bam}' > ${meta.id}.stats
    samtools flagstat --input-fmt-option reference='${genome_fasta}' -@ ${task.cpus} '${bam}' > ${meta.id}.flagstat
    samtools idxstats '${bam}' > ${meta.id}.idxstats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.stats ${meta.id}.flagstat ${meta.id}.idxstats
    echo '"${task.process}": {samtools: stub}' > versions.yml
    """
}

process SAMTOOLS_CRAM_TO_BAM {
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(cram)
    path  genome_fasta

    output:
    tuple val(meta), path('*.bam'), emit: bam
    path 'versions.yml'           , emit: versions

    script:
    """
    samtools view -@ ${task.cpus} -b -T '${genome_fasta}' -o ${meta.id}.bam '${cram}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.bam
    echo '"${task.process}": {samtools: stub}' > versions.yml
    """
}

process INFER_STRANDEDNESS {
    label 'process_medium'
    tag "${meta.id}"

    // Deliberately takes an UNINDEXED alignment. Strandedness has to be
    // known before the real alignment (it sets --rna-strandness, which drives
    // the XS tag StringTie depends on), so this runs on a cheap unstranded
    // alignment of a read subsample and streams it with fetch(until_eof).
    // The probe is CRAM, so pysam must be handed the reference to decode it.
    input:
    tuple val(meta), path(bam)
    tuple path(fasta), path(fai)

    output:
    tuple val(meta), path('*.strandedness.json'), emit: json
    path 'versions.yml'                         , emit: versions

    script:
    def strict = params.strict_strandedness ? '--strict' : ''
    """
    infer_strandedness.py \\
        --bam '${bam}' --fasta '${fasta}' \\
        --sample-id '${meta.id}' \\
        --declared '${meta.declared ?: meta.strandedness}' \\
        --min-reads ${params.strandedness_min_reads} \\
        --threshold ${params.strandedness_threshold} \\
        ${strict} \\
        --out ${meta.id}.strandedness.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pysam: \$(python3 -c 'import pysam; print(pysam.__version__)')
    END_VERSIONS
    """

    stub:
    """
    echo '{"sample_id":"${meta.id}","inferred":"unstranded","effective":"unstranded","conflict":false,"fraction_agreeing":0.5,"informative_reads":0}' > ${meta.id}.strandedness.json
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
