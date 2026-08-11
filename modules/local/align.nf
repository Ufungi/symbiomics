//
// Alignment: HISAT2 -> coordinate-sorted BAM -> index -> stats.
//
// Two things here are load-bearing and were bugs in the pipeline this replaces:
//   * SAM is never written to disk; hisat2 pipes straight into samtools sort.
//   * The index is CSI whenever any contig exceeds the BAI 512 Mb ceiling.
//     BAI cannot address a 2 Gb conifer chromosome at all.
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

    // hisat2's wrapper decompresses gzipped reads through FIFOs at a hardcoded
    // /tmp/$$.inpipe* (it ignores TMPDIR, and its --temp-directory flag is
    // forwarded to the aligner binary, which rejects it). Singularity's --pid
    // namespace gives every container a small $$, so concurrent tasks -- and
    // reruns after a crash left a stale FIFO behind -- collide on the same path
    // in the shared host /tmp. Giving each task its own /tmp removes the class
    // of bug rather than the instance.
    containerOptions = {
        workflow.containerEngine in ['singularity','apptainer'] ? "--bind ${task.workDir}:/tmp"
      : workflow.containerEngine == 'docker'                    ? "-v ${task.workDir}:/tmp"
      : ''
    }

    input:
    tuple val(meta), path(reads)
    path  index
    val   index_prefix
    val   use_mmap

    output:
    tuple val(meta), path('*.unsorted.bam'), emit: bam
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

    def input_args
    if( meta.single_end ) {
        input_args = "-U ${reads.join(',')}"
    } else {
        def half = (reads.size() / 2) as int
        input_args = "-1 ${reads[0..<half].join(',')} -2 ${reads[half..<reads.size()].join(',')}"
    }
    """
    hisat2 -p ${task.cpus} -x ${index}/${index_prefix} ${input_args} \\
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
    touch ${meta.id}.unsorted.bam ${meta.id}.summary.txt
    echo '"${task.process}": {hisat2: stub}' > versions.yml
    """
}

process SAMTOOLS_SORT {
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*.sorted.bam'), emit: bam
    path 'versions.yml'                  , emit: versions

    script:
    """
    samtools sort -@ ${task.cpus} -m ${params.sort_mem} \\
        -T ${meta.id}.tmp -o ${meta.id}.sorted.bam '${bam}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sorted.bam
    echo '"${task.process}": {samtools: stub}' > versions.yml
    """
}

process SAMTOOLS_INDEX {
    label 'process_low'
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam)
    val   index_type          // csi | bai | both

    output:
    tuple val(meta), path(bam), path('*.{csi,bai}'), emit: bam
    path 'versions.yml'                            , emit: versions

    script:
    // CSI is always safe; BAI is only added when every contig is addressable,
    // because some downstream tools still hardcode a .bai suffix.
    def make_csi = index_type in ['csi', 'both']
    def make_bai = index_type in ['bai', 'both']
    """
    ${make_csi ? "samtools index -c -@ ${task.cpus} '${bam}'" : ''}
    ${make_bai ? "samtools index -b -@ ${task.cpus} '${bam}'" : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${bam}.csi
    echo '"${task.process}": {samtools: stub}' > versions.yml
    """
}

process BAM_STATS {
    label 'process_low'
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(index)

    output:
    path '*.stats'      , emit: stats
    path '*.flagstat'   , emit: flagstat
    path '*.idxstats'   , emit: idxstats
    path 'versions.yml' , emit: versions

    script:
    """
    samtools stats    -@ ${task.cpus} '${bam}' > ${meta.id}.stats
    samtools flagstat -@ ${task.cpus} '${bam}' > ${meta.id}.flagstat
    samtools idxstats '${bam}'                 > ${meta.id}.idxstats

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

process INFER_STRANDEDNESS {
    label 'process_medium'
    tag "${meta.id}"

    // Deliberately takes an UNSORTED, UNINDEXED bam. Strandedness has to be
    // known before the real alignment (it sets --rna-strandness, which drives
    // the XS tag StringTie depends on), so this runs on a cheap unstranded
    // alignment of a read subsample and streams it with fetch(until_eof).
    // Sorting it first would be pure waste.
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
