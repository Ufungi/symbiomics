//
// Parse and validate the samplesheet, then build the per-sample channel.
//
// Validation runs as a real process (containerised, versioned, cached) rather
// than in Groovy, so the resolved samplesheet is a first-class, reproducible
// artefact that can be fed straight back in as --input.
//

process INPUT_CHECK {
    label 'process_single'
    tag "${samplesheet.name}"

    input:
    path samplesheet
    val  raw_dir
    val  default_strandedness
    val  default_genome_id
    val  only_sample

    output:
    path 'samplesheet.resolved.tsv' , emit: tsv
    path 'samplesheet.resolved.json', emit: json
    path 'versions.yml'             , emit: versions

    script:
    def raw     = raw_dir     ? "--raw-dir '${raw_dir}'"          : ''
    def gid     = default_genome_id ? "--default-genome-id '${default_genome_id}'" : ''
    def only    = only_sample ? "--sample '${only_sample}'"       : ''
    """
    infer_layout.py \\
        --input '${samplesheet}' \\
        --out-tsv samplesheet.resolved.tsv \\
        --out-json samplesheet.resolved.json \\
        --default-strandedness '${default_strandedness}' \\
        ${raw} ${gid} ${only}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    // Run the real parser but skip filesystem checks. A stub that invented its
    // own fake sample would validate nothing; this way -stub-run still exercises
    // the header handling and the SE/PE inference ladder.
    """
    infer_layout.py \\
        --input '${samplesheet}' \\
        --out-tsv samplesheet.resolved.tsv \\
        --out-json samplesheet.resolved.json \\
        --default-strandedness '${default_strandedness}' \\
        --no-check-files ${raw_dir ? "--raw-dir '${raw_dir}'" : ''}
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}

workflow INPUT_CHECK_WF {
    take:
    samplesheet          // path
    raw_dir              // val, may be null
    default_strandedness // val
    default_genome_id    // val, may be null
    only_sample          // val, may be null

    main:
    INPUT_CHECK(samplesheet, raw_dir ?: '', default_strandedness,
                default_genome_id ?: '', only_sample ?: '')

    reads = INPUT_CHECK.out.tsv
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = [
                id          : row.sample_id,
                single_end  : row.layout == 'SE',
                layout      : row.layout,
                strandedness: row.strandedness,
                condition   : row.condition,
                replicate   : row.replicate,
                read_type   : row.read_type,
                genome_id   : row.genome_id == '-' ? null : row.genome_id,
                use_for     : row.use_for.tokenize(','),
                prealigned  : row.bam != '-',
            ]
            // infer_layout.py has already verified existence; re-checking here
            // would only duplicate that, and must not fire during -stub-run.
            def must_exist = !workflow.stubRun
            if( meta.prealigned ) {
                return [ meta, [], file(row.bam, checkIfExists: must_exist) ]
            }
            def r1 = row.fastq_1.tokenize(',').collect { file(it, checkIfExists: must_exist) }
            def r2 = row.fastq_2 == '-' ? [] : row.fastq_2.tokenize(',').collect { file(it, checkIfExists: must_exist) }
            return [ meta, r1 + r2, null ]
        }

    // Split pre-aligned rows out: they skip trim + align entirely.
    reads
        .branch { meta, fastqs, bam ->
            prealigned: meta.prealigned
                return [ meta, bam ]
            to_align: true
                return [ meta, fastqs ]
        }
        .set { split }

    emit:
    reads      = split.to_align       // [ meta, [fastq...] ]
    prealigned = split.prealigned     // [ meta, bam ]
    resolved   = INPUT_CHECK.out.tsv
    versions   = INPUT_CHECK.out.versions
}
