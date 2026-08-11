//
// Read-level processes: QC, subsampling, trimming.
//
// Channel contract for everything in this file:
//     input  [ meta, [ fastq... ] ]   1 file for SE, 2 for PE (already paired
//                                     and, for technical replicates, in
//                                     R1a,R1b,R2a,R2b order)
//     output [ meta, [ fastq... ] ]   same shape
//

process FASTQC {
    label 'process_low'
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.zip') , emit: zip
    path '*.html'                  , emit: html
    path 'versions.yml'            , emit: versions

    script:
    // FastQC names its output after the input file. Two samples that point at
    // the same fastq -- a technical replicate, or a re-run under a second
    // sample_id -- would then emit identically named reports and collide when
    // MultiQC stages them together. Prefix everything with the sample id.
    """
    fastqc --threads ${task.cpus} --quiet ${reads}

    for f in *_fastqc.zip *_fastqc.html; do
        [ -e "\$f" ] || continue
        case "\$f" in
            ${meta.id}_*) ;;
            *) mv "\$f" "${meta.id}_\$f" ;;
        esac
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_stub_fastqc.zip ${meta.id}_stub_fastqc.html
    echo '"${task.process}": {fastqc: stub}' > versions.yml
    """
}

process SEQKIT_SUBSAMPLE {
    label 'process_low'
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)
    val  n_reads
    val  sampler
    val  seed

    output:
    tuple val(meta), path('sub_*.fastq.gz'), emit: reads
    path 'versions.yml'                    , emit: versions

    script:
    // `head` is pair-safe by construction (same first N records on both mates)
    // and needs no seed; `sample` reproduces the lab's previous seqtk-style
    // behaviour but requires identical record order in R1/R2.
    def cmd = sampler == 'sample' ? "seqkit sample -s ${seed} -n ${n_reads}"
                                  : "seqkit head -n ${n_reads}"
    """
    i=1
    for fq in ${reads}; do
        ${cmd} "\$fq" | gzip -c > "sub_\${i}_${meta.id}.fastq.gz"
        i=\$((i+1))
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit v//')
    END_VERSIONS
    """

    stub:
    """
    i=1; for fq in ${reads}; do echo | gzip -c > "sub_\${i}_${meta.id}.fastq.gz"; i=\$((i+1)); done
    echo '"${task.process}": {seqkit: stub}' > versions.yml
    """
}

process FASTP {
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('trim_*.fastq.gz'), emit: reads
    tuple val(meta), path('*.fastp.json')   , emit: json
    path '*.fastp.html'                     , emit: html
    path 'versions.yml'                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    if( meta.single_end ) {
        """
        cat ${reads} > merged_1.fq.gz
        fastp --in1 merged_1.fq.gz --out1 trim_1_${meta.id}.fastq.gz \\
            --thread ${task.cpus} ${args} \\
            --json ${meta.id}.fastp.json --html ${meta.id}.fastp.html
        rm -f merged_1.fq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
        END_VERSIONS
        """
    } else {
        // reads arrive as [R1a, R1b, ..., R2a, R2b, ...]; concatenating each
        // half preserves record order across technical replicates.
        """
        n=\$(( \$(ls -1 ${reads} | wc -l) / 2 ))
        ls -1 ${reads} | head -n \$n | xargs cat > merged_1.fq.gz
        ls -1 ${reads} | tail -n \$n | xargs cat > merged_2.fq.gz

        fastp --in1 merged_1.fq.gz --in2 merged_2.fq.gz \\
            --out1 trim_1_${meta.id}.fastq.gz --out2 trim_2_${meta.id}.fastq.gz \\
            --detect_adapter_for_pe --thread ${task.cpus} ${args} \\
            --json ${meta.id}.fastp.json --html ${meta.id}.fastp.html
        rm -f merged_1.fq.gz merged_2.fq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
        END_VERSIONS
        """
    }

    stub:
    """
    echo | gzip -c > trim_1_${meta.id}.fastq.gz
    ${meta.single_end ? '' : "echo | gzip -c > trim_2_${meta.id}.fastq.gz"}
    echo '{}' > ${meta.id}.fastp.json
    touch ${meta.id}.fastp.html
    echo '"${task.process}": {fastp: stub}' > versions.yml
    """
}
