//
// Genome statistics and the policy engine.
//
// GENOME_STATS streams the FASTA once; DECIDE_STRATEGY turns those numbers plus
// the available evidence and the user's parameters into strategy.yml, which the
// rest of the pipeline reads instead of re-deriving decisions in five places.
//

process SAMTOOLS_FAIDX {
    label 'process_single'
    tag "${fasta.name}"

    input:
    path fasta

    output:
    path "${fasta}.fai" , emit: fai
    path 'versions.yml' , emit: versions

    script:
    """
    samtools faidx '${fasta}' --fai-idx '${fasta}.fai'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    printf 'chr1\\t1000\\t6\\t60\\t61\\n' > ${fasta}.fai
    echo '"${task.process}": {samtools: stub}' > versions.yml
    """
}

process GENOME_STATS {
    label 'process_low'
    tag "${meta.id}"

    input:
    tuple val(meta), path(fasta)
    val   sample_bp

    output:
    tuple val(meta), path('genome_stats.json'), emit: stats
    path 'versions.yml'                       , emit: versions

    script:
    """
    genome_stats.py --fasta '${fasta}' --out genome_stats.json --sample-bp ${sample_bp}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    cat > genome_stats.json <<'JSON'
    {"total_bp": 100000, "contigs": 1, "max_contig_bp": 100000,
     "n_fraction": 0.0, "softmask_fraction": 0.5,
     "composition_scanned_bp": 100000, "composition_is_sampled": false}
    JSON
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}

process DECIDE_STRATEGY {
    label 'process_single'
    tag "${meta.id}"

    input:
    tuple val(meta), path(stats)
    val  evidence      // map: has_rna / has_protein / has_isoseq / n_samples
    val  index_bytes

    output:
    tuple val(meta), path('strategy.json'), emit: json
    path 'strategy.yml'                   , emit: yml
    path 'versions.yml'                   , emit: versions

    script:
    def species  = meta.species     ? "--species '${meta.species}'"       : ''
    def clade    = params.clade     ? "--clade '${params.clade}'"         : ''
    def index    = params.hisat2_index ? "--hisat2-index '${params.hisat2_index}'" : ''
    def repeatlib = params.repeat_lib  ? "--repeat-lib '${params.repeat_lib}'"     : ''
    def budget   = params.time_budget  ? "--time-budget '${params.time_budget}'"   : ''
    """
    decide_strategy.py \\
        --genome-stats '${stats}' \\
        --genome-id '${meta.id}' ${species} ${clade} \\
        --taxon '${params.taxon}' \\
        --size-class '${params.genome_size_class}' \\
        --threshold-tiny ${params.size_threshold_tiny} \\
        --threshold-small ${params.size_threshold_small} \\
        --threshold-medium ${params.size_threshold_medium} \\
        --threshold-large ${params.size_threshold_large} \\
        --masker '${params.masker}' \\
        --premasked ${params.premasked} ${repeatlib} \\
        --aligner '${params.aligner}' ${index} \\
        --hisat2-index-bytes ${index_bytes} \\
        --bam-index '${params.bam_index}' \\
        --has-rna ${evidence.has_rna} \\
        --has-protein ${evidence.has_protein} \\
        --has-isoseq ${evidence.has_isoseq} \\
        --n-samples ${evidence.n_samples} \\
        --priority '${params.priority}' ${budget} \\
        --gpu ${params.gpu == 'auto' ? 'false' : params.gpu} \\
        --max-cpus ${params.max_cpus} \\
        --out strategy.yml --out-json strategy.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo '{"genome":{"id":"${meta.id}","size_class":"tiny","taxon":"other"},"decisions":{},"budget":{},"notes":[]}' > strategy.json
    echo 'genome: {id: ${meta.id}}' > strategy.yml
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}

workflow GENOME_PREP {
    take:
    genome        // [ meta, fasta ]
    evidence      // val map
    index_bytes   // val long

    main:
    ch_versions = Channel.empty()

    SAMTOOLS_FAIDX(genome.map { _meta, fasta -> fasta })
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    GENOME_STATS(genome, params.test_mode ? 20_000_000 : 200_000_000)
    ch_versions = ch_versions.mix(GENOME_STATS.out.versions)

    DECIDE_STRATEGY(GENOME_STATS.out.stats, evidence, index_bytes)
    ch_versions = ch_versions.mix(DECIDE_STRATEGY.out.versions)

    // Decode strategy.json once and hand the map to every downstream consumer.
    //
    // JsonSlurperClassic, not JsonSlurper: the latter returns a LazyMap that
    // builds its backing store on first access, which races and NPEs when the
    // same map is read from several GPars dataflow actors at once.
    // JsonSlurperClassic returns plain LinkedHashMaps.
    strategy = DECIDE_STRATEGY.out.json
        .map { meta, json -> [ meta, new groovy.json.JsonSlurperClassic().parseText(json.text) ] }

    emit:
    fai      = SAMTOOLS_FAIDX.out.fai
    stats    = GENOME_STATS.out.stats
    strategy = strategy                   // [ meta, Map ]
    yml      = DECIDE_STRATEGY.out.yml
    versions = ch_versions
}
