//
// Transcript assembly and counting.
//
// Both counters carry a fix for a real bug in the pipeline this replaces:
//   * HTSeq is always given `-r pos`. Its default is `-r name`, and running it
//     on coordinate-sorted BAM silently loses PE mates and balloons memory.
//   * Counting is over `exon` features grouped by gene_id, never over `gene`.
//     Counting the whole gene span pulls in intronic reads.
//

process STRINGTIE_ASSEMBLE {
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(index)
    path  guide_gff

    output:
    tuple val(meta), path('*.transcripts.gtf'), emit: gtf
    path 'versions.yml'                       , emit: versions

    script:
    def strand = meta.strandedness == 'forward' ? '--fr'
               : meta.strandedness == 'reverse' ? '--rf' : ''
    def guide  = guide_gff.name != 'NO_FILE' ? "-G ${guide_gff}" : ''
    def extra  = params.stringtie_extra ?: ''
    """
    stringtie '${bam}' -p ${task.cpus} ${strand} ${guide} ${extra} \\
        -o ${meta.id}.transcripts.gtf -l ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stringtie: \$(stringtie --version)
    END_VERSIONS
    """

    stub:
    """
    echo '# stringtie stub' > ${meta.id}.transcripts.gtf
    echo '"${task.process}": {stringtie: stub}' > versions.yml
    """
}

process STRINGTIE_MERGE {
    label 'process_medium'
    tag 'merge'

    input:
    path gtfs
    path guide_gff

    output:
    path 'merged.gtf'   , emit: gtf
    path 'versions.yml' , emit: versions

    script:
    // stringtie --merge aborts on the first input with no transcripts, which
    // would let one low-coverage sample kill the merge for everything else.
    // Drop empty GTFs, say which, and emit an empty merged.gtf if none survive
    // rather than failing -- the emptiness is then visible downstream instead
    // of appearing as a crash here.
    def guide = guide_gff.name != 'NO_FILE' ? "-G ${guide_gff}" : ''
    """
    # awk, not grep -E: GNU grep does not read \\t as a tab in a BRE/ERE, so a
    # grep-based test silently matches nothing and would discard every input.
    : > mergelist.txt
    for gtf in ${gtfs}; do
        if awk -F'\\t' '\$1 !~ /^#/ && \$3 == "transcript" { found = 1; exit } END { exit !found }' "\$gtf"; then
            echo "\$gtf" >> mergelist.txt
        else
            echo "WARN  ~ [stringtie] \$gtf has no transcripts; excluded from the merge" >&2
        fi
    done

    if [ -s mergelist.txt ]; then
        stringtie --merge -p ${task.cpus} ${guide} -o merged.gtf mergelist.txt
    else
        echo "WARN  ~ [stringtie] no sample yielded transcripts; merged.gtf is empty" >&2
        printf '# stringtie merge: no input GTF contained transcripts\\n' > merged.gtf
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stringtie: \$(stringtie --version)
    END_VERSIONS
    """

    stub:
    """
    echo '# merged stub' > merged.gtf
    echo '"${task.process}": {stringtie: stub}' > versions.yml
    """
}

process INSPECT_ANNOTATION {
    label 'process_single'
    tag "${annotation.name}"

    input:
    path annotation

    output:
    path 'annotation_spec.json', emit: spec
    path 'versions.yml'        , emit: versions

    script:
    """
    inspect_annotation.py \\
        --annotation '${annotation}' \\
        --want-feature '${params.count_feature}' \\
        --want-attribute '${params.count_attribute}' \\
        --out annotation_spec.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo '{"feature":"exon","attribute":"gene_id","n_groups":0,"notes":[]}' > annotation_spec.json
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}

process FEATURECOUNTS {
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(index)
    path  annotation
    val   spec          // [feature: ..., attribute: ...] from INSPECT_ANNOTATION

    output:
    tuple val(meta), path('*.featureCounts.txt')        , emit: counts
    tuple val(meta), path('*.featureCounts.txt.summary'), emit: summary
    path 'versions.yml'                                 , emit: versions

    script:
    def s = meta.strandedness == 'forward' ? 1 : meta.strandedness == 'reverse' ? 2 : 0
    def paired = meta.single_end ? '' : '-p --countReadPairs'
    def multi  = params.count_multimappers ? '-M' : ''
    """
    # An annotation with no countable features makes featureCounts exit 255 with
    # a wall of banner text. Say what is actually wrong instead.
    if ! awk -F'\\t' -v t=${spec.feature} \\
         '\$1 !~ /^#/ && \$3 == t { found = 1; exit } END { exit !found }' '${annotation}'; then
        echo "ERROR ~ [quantify] ${annotation} contains no '${spec.feature}' features." >&2
        echo "        Nothing can be counted against it. If this is a StringTie merge," >&2
        echo "        no sample yielded transcripts; check 45_assemble/ and the alignment rates." >&2
        exit 1
    fi

    featureCounts \\
        -T ${task.cpus} ${paired} ${multi} \\
        -t ${spec.feature} -g ${spec.attribute} \\
        -s ${s} -F GTF \\
        -a '${annotation}' -o ${meta.id}.featureCounts.txt '${bam}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$(featureCounts -v 2>&1 | grep -oE 'v[0-9.]+' | tr -d 'v')
    END_VERSIONS
    """

    stub:
    """
    printf 'Geneid\\tChr\\tStart\\tEnd\\tStrand\\tLength\\t${meta.id}\\ng1\\tc\\t1\\t2\\t+\\t2\\t0\\n' > ${meta.id}.featureCounts.txt
    touch ${meta.id}.featureCounts.txt.summary
    echo '"${task.process}": {subread: stub}' > versions.yml
    """
}

process HTSEQ_COUNT {
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(index)
    path  annotation
    val   spec

    output:
    tuple val(meta), path('*.htseq.tsv'), emit: counts
    path 'versions.yml'                 , emit: versions

    script:
    def s = meta.strandedness == 'forward' ? 'yes'
          : meta.strandedness == 'reverse' ? 'reverse' : 'no'
    """
    # -r pos is mandatory here: the BAM is coordinate-sorted and HTSeq defaults
    # to -r name, which drops PE mates and blows up memory.
    htseq-count \\
        --order pos --stranded ${s} \\
        --type ${spec.feature} --idattr ${spec.attribute} \\
        --nprocesses ${task.cpus} \\
        --counts_output ${meta.id}.htseq.tsv \\
        '${bam}' '${annotation}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htseq: \$(htseq-count --version)
    END_VERSIONS
    """

    stub:
    """
    printf 'g1\\t0\\n' > ${meta.id}.htseq.tsv
    echo '"${task.process}": {htseq: stub}' > versions.yml
    """
}

process MERGE_COUNTS {
    label 'process_low'
    tag "${tool}"

    input:
    path  count_files
    val   tool          // featurecounts | htseq
    path  strandedness_json

    output:
    path "counts_${tool}.tsv"  , emit: matrix
    path "counts_${tool}.stats.tsv", emit: stats
    path 'versions.yml'        , emit: versions

    script:
    def strand_arg = strandedness_json.name != 'NO_FILE' ? "--strandedness ${strandedness_json}" : ''
    """
    merge_counts.py --tool ${tool} ${strand_arg} \\
        --out counts_${tool}.tsv --stats counts_${tool}.stats.tsv \\
        ${count_files}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf 'gene_id\\n' > counts_${tool}.tsv
    printf 'metric\\tvalue\\n' > counts_${tool}.stats.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}

process MULTIQC {
    label 'process_low'
    tag 'report'

    input:
    path  files, stageAs: 'input/*'
    path  config

    output:
    path 'multiqc_report.html' , emit: report
    path 'multiqc_report_data' , emit: data
    path 'versions.yml'        , emit: versions

    script:
    // --filename is explicit because MultiQC otherwise derives the output name
    // from `title` in the config, which would drift away from what this process
    // declares as its output.
    def cfg = config.name != 'NO_FILE' ? "--config ${config}" : ''
    """
    multiqc --force ${cfg} --filename multiqc_report.html input/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //')
    END_VERSIONS
    """

    stub:
    """
    touch multiqc_report.html && mkdir -p multiqc_report_data
    echo '"${task.process}": {multiqc: stub}' > versions.yml
    """
}
