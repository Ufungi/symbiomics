//
// Multi-genome read mapping: the same trimmed reads aligned against N
// candidate reference genomes (e.g. host + one or more symbiont genomes) so
// per-genome mapping rates can be compared directly.
//
// Deliberately separate from modules/local/align.nf rather than reusing
// HISAT2_BUILD/HISAT2_ALIGN by alias: those processes are broadcast-index
// shaped (one index, many samples) because the main arm only ever aligns
// against a single genome. Here the index itself varies per task, so it has
// to travel through the tuple instead of arriving as a scalar `path`/`val`
// input -- that changes both processes' input shape, not just their names.
//

process HISAT2_BUILD_MULTI {
    label 'process_high'
    tag "${genome_meta.id}"

    input:
    tuple val(genome_meta), path(fasta)

    output:
    tuple val(genome_meta), path('hisat2_index'), emit: index
    path 'versions.yml'                         , emit: versions

    script:
    def large = genome_meta.large_index ? '--large-index' : ''
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

process HISAT2_ALIGN_MULTI {
    label 'process_high'
    tag "${meta.id}"

    // Same /tmp-collision fix as HISAT2_ALIGN -- see modules/local/align.nf.
    containerOptions = {
        workflow.containerEngine in ['singularity','apptainer'] ? "--bind ${task.workDir}:/tmp"
      : workflow.containerEngine == 'docker'                    ? "-v ${task.workDir}:/tmp"
      : ''
    }

    input:
    tuple val(meta), path(reads), path(index)

    output:
    tuple val(meta), path('*.unsorted.bam'), emit: bam
    tuple val(meta), path('*.summary.txt') , emit: summary
    path 'versions.yml'                    , emit: versions

    script:
    // mmap follows the genome this task aligns against, not a single global
    // choice -- a run mixing one huge candidate genome with several small
    // ones should only pay the shared-index-cache cost where it matters.
    def mmap  = meta.large_index ? '--mm' : ''
    def extra = params.hisat2_extra ?: ''
    def strand = ''
    if( meta.strandedness == 'forward' ) strand = meta.single_end ? '--rna-strandness F' : '--rna-strandness FR'
    if( meta.strandedness == 'reverse' ) strand = meta.single_end ? '--rna-strandness R' : '--rna-strandness RF'

    def input_args
    if( meta.single_end ) {
        input_args = "-U ${reads.join(',')}"
    } else {
        def half = (reads.size() / 2) as int
        input_args = "-1 ${reads[0..<half].join(',')} -2 ${reads[half..<reads.size()].join(',')}"
    }
    """
    hisat2 -p ${task.cpus} -x ${index}/genome ${input_args} \\
        ${strand} ${mmap} ${extra} \\
        --new-summary --summary-file ${meta.id}.summary.txt \\
        | samtools view -@ 2 -bS -o ${meta.id}.unsorted.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$(hisat2 --version | head -1 | sed 's/.*version //')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.unsorted.bam
    printf 'Overall alignment rate: 0.00%%\\n' > ${meta.id}.summary.txt
    echo '"${task.process}": {hisat2: stub}' > versions.yml
    """
}

process SUMMARIZE_MULTI_GENOME_MAPPING {
    label 'process_single'

    input:
    path summaries  // *.summary.txt, one per (sample, genome)
    path manifest   // TSV: sample_id, genome_id, summary_file

    output:
    path 'mapping_matrix.tsv'  , emit: matrix
    path 'best_genome.tsv'     , emit: best_genome
    path 'mapping_summary.json', emit: json
    path 'versions.yml'        , emit: versions

    script:
    """
    summarize_multi_genome_mapping.py \\
        --manifest '${manifest}' \\
        --ambiguous-margin ${params.mapping_ambiguous_margin} \\
        --out-matrix mapping_matrix.tsv \\
        --out-best best_genome.tsv \\
        --out-json mapping_summary.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch mapping_matrix.tsv best_genome.tsv mapping_summary.json
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
