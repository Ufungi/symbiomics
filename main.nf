#!/usr/bin/env nextflow
/*
 * symbiomics -- decision-driven eukaryote genome annotation + mRNA-seq pipeline
 * https://github.com/Ufungi/symbiomics
 *
 * Entry points:
 *   (default)          full pipeline
 *   -entry strategy    resolve and print the execution plan, run nothing else
 *   -entry preflight   environment checks only
 */

nextflow.enable.dsl = 2

include { SYMBIOMICS     } from './workflows/symbiomics'
include { FUNCTIONAL   } from './workflows/functional'
include { GENOME_PREP  } from './subworkflows/local/genome_prep'
include { INPUT_CHECK_WF } from './subworkflows/local/input_check'
include { HAPLOTYPE_PAIRING } from './subworkflows/local/haplotype_pairing'
include { INSPECT_VARIATION_FILE } from './modules/local/variation'
include { MULTI_GENOME_MAPPING } from './subworkflows/local/multi_genome_mapping'
include { MULTIQC } from './modules/local/quantify'

def helpMessage() {
    log.info """
    ${workflow.manifest.name} ${workflow.manifest.version}
    ${workflow.manifest.description}

    Usage:
      scripts/symbiomics run . -profile singularity,local64 \\
          --input samplesheet.tsv --genome genome.fasta --outdir results

    Required:
      --input                RNA-seq samplesheet (TSV or CSV). See assets/samplesheet.example.tsv
      --genome               genome FASTA (alias: --fasta)

    Commonly used:
      --outdir               output directory                    [${params.outdir}]
      --raw_dir              base directory for relative fastq paths
      --taxon                auto|fungi|plant|animal|other        [${params.taxon}]
      --genome_size_class    auto|tiny|small|medium|large|huge    [${params.genome_size_class}]
      --strandedness         auto|unstranded|forward|reverse      [${params.strandedness}]
      --hisat2_index         reuse a prebuilt HISAT2 index prefix
      --premasked            genome is already softmasked         [${params.premasked}]
      --counter              featurecounts|htseq|both|none        [${params.counter}]
      --steps / --skip_steps comma lists of stage names
      --sample               run a single sample_id
      --test_mode            subsample inputs                     [${params.test_mode}]

    Planning:
      --dry_run_strategy     resolve the execution plan and stop. Do this before
                             launching anything that will run for days.

    Multi-genome symbiont mapping (-entry symbiont_mapping):
      --input                RNA-seq samplesheet (as above)
      --mapping_genomes      TSV: genome_id, fasta, taxon -- 2+ candidate
                             genomes (e.g. host + symbiont) to score the same
                             reads against. See assets/mapping_genomes.example.tsv

    Profiles: singularity, docker, conda, gpu, local64, test, test_fungus, test_mapping
    """.stripIndent()
}

workflow {
    if( params.help ) {
        helpMessage()
        return
    }
    // --genome is not required for the two paths that operate on sequences
    // directly: --transcript_fasta (Salmon-only quantification, no alignment)
    // and --proteome (functional annotation, no structural annotation at all).
    def genome_free_ok = params.transcript_fasta || params.proteome
    if( !(params.genome ?: params.fasta) && !genome_free_ok ) {
        helpMessage()
        error "Missing required parameter: --genome (or --fasta) -- " +
              "unless running --quant_engine salmon with --transcript_fasta, " +
              "or -entry functional with --proteome."
    }
    if( !params.input && !params.dry_run_strategy && !params.proteome ) {
        helpMessage()
        error "Missing required parameter: --input (a samplesheet). " +
              "Use --dry_run_strategy to resolve a plan without one, or " +
              "-entry functional --proteome for protein-only functional annotation."
    }
    SYMBIOMICS()
}

/*
 * Resolve the execution plan without running the pipeline. This exists because
 * a huge-genome run is a multi-week commitment and the plan should be reviewed
 * before it starts, not discovered from the log afterwards.
 */
workflow strategy {
    genome_file = file(params.genome ?: params.fasta, checkIfExists: true)
    genome_id   = params.genome_id ?: genome_file.simpleName
    genome_meta = [ id: genome_id, species: params.species ?: genome_id ]

    if( params.input ) {
        INPUT_CHECK_WF(file(params.input, checkIfExists: true), params.raw_dir,
                       params.strandedness, genome_id, params.sample)
        evidence = INPUT_CHECK_WF.out.resolved.map { tsv ->
            def rows = tsv.readLines().findAll { it && !it.startsWith('sample_id') }
            [ has_rna: rows.size() > 0,
              has_protein: (params.proteins ?: '') != '' && params.proteins != 'auto',
              has_isoseq: rows.any { it.split('\t').size() > 7 && it.split('\t')[7] != 'short' },
              n_samples: rows.size() ]
        }
    } else {
        evidence = Channel.value([ has_rna: false, has_protein: false,
                                   has_isoseq: false, n_samples: 0 ])
    }

    index_bytes = params.hisat2_index
        ? file(params.hisat2_index).parent.listFiles()
              .findAll { it.name.endsWith('.ht2') || it.name.endsWith('.ht2l') }
              .sum { it.size() } ?: 0
        : 0

    GENOME_PREP(Channel.value([ genome_meta, genome_file ]), evidence, index_bytes)
    GENOME_PREP.out.yml.view { it.text }
}

/*
 * -entry functional -- protein-in functional annotation. No genome, no
 * samplesheet: scripts/symbiomics run . -entry functional --proteome p.faa
 */
workflow functional {
    FUNCTIONAL()
}

/*
 * -entry pairing -- HA<->HB assembly alignment (minimap2 + SyRI): a phased
 * VCF for phASER and a gene-level allele-pairing table for the diploid Salmon
 * EM split. Standalone so this expensive step runs once per genome pair and
 * both downstream consumers (F6, F7) reuse its output.
 */
workflow pairing {
    if( !params.genome_hb ) {
        error "Missing required parameter: --genome_hb <second-haplotype genome fasta>"
    }
    ha = file(params.genome ?: params.fasta, checkIfExists: true)
    hb = file(params.genome_hb, checkIfExists: true)
    // Distinct sentinels: both can be "not given" at once, and staging two
    // files literally named the same thing into one task dir collides.
    gff_ha = params.gff3_ha ? file(params.gff3_ha, checkIfExists: true) : file("${projectDir}/assets/NO_FILE_GFF3_HA")
    gff_hb = params.gff3_hb ? file(params.gff3_hb, checkIfExists: true) : file("${projectDir}/assets/NO_FILE_GFF3_HB")

    HAPLOTYPE_PAIRING(
        Channel.value(ha), Channel.value(hb),
        Channel.value(gff_ha), Channel.value(gff_hb),
        params.minimap_preset
    )
    HAPLOTYPE_PAIRING.out.phased_vcf.view   { "phased VCF: ${it}" }
    HAPLOTYPE_PAIRING.out.allele_table.view { "allele pairing table: ${it}" }
}

/*
 * -entry variation -- route 2 for the ASE phased VCF: probe the dataset's own
 * resequencing variant file rather than compute one from HA/HB alignment.
 * See docs/allele_specific_expression.md.
 */
workflow variation {
    if( !params.resequencing_variation_file ) {
        error "Missing required parameter: --resequencing_variation_file <path>"
    }
    INSPECT_VARIATION_FILE(
        file(params.resequencing_variation_file, checkIfExists: true),
        params.resequencing_accession ?: ''
    )
    INSPECT_VARIATION_FILE.out.vcf.view { "VCF: ${it}" }
}

/*
 * -entry symbiont_mapping -- map ONE set of transcriptome reads against
 * SEVERAL candidate genomes (e.g. host + one or more symbiont genomes, or
 * several candidate symbiont species of unknown identity) and compare
 * per-genome mapping rates. This is what lets a mixed or ambiguously-sourced
 * RNA-seq sample be assigned to the genome it actually came from, rather
 * than committing to one genome before running anything. See
 * docs/symbiont_mapping.md.
 */
workflow symbiont_mapping {
    if( !params.input ) {
        error "Missing required parameter: --input <RNA-seq samplesheet>"
    }
    if( !params.mapping_genomes ) {
        error "Missing required parameter: --mapping_genomes <tsv: genome_id, fasta, taxon> " +
              "-- see assets/mapping_genomes.example.tsv"
    }

    genomes_file = file(params.mapping_genomes, checkIfExists: true)
    rows = genomes_file.readLines()
        .findAll { it && !it.startsWith('#') && !it.startsWith('genome_id') }
    if( rows.size() < 2 ) {
        error "--mapping_genomes needs at least 2 genomes to compare against " +
              "(found ${rows.size()} in ${genomes_file}) -- a single genome has " +
              "nothing to be compared to; use the default entry point instead."
    }

    // Parsed synchronously, not inside a channel `.map` -- the table is a
    // param known before the pipeline starts (like --genome_hb in `pairing`),
    // so there is nothing gained by deferring it into the dataflow graph, and
    // duplicate-id tracking is far simpler as a plain Groovy loop than as
    // mutable state captured by a reactive closure.
    def seen_ids = [] as Set
    def genome_rows = rows.collect { line ->
        def f = line.split('\t')
        if( f.size() < 2 ) {
            error "--mapping_genomes: malformed row (need at least genome_id, fasta): '${line}'"
        }
        def gid = f[0]
        if( !seen_ids.add(gid) ) {
            error "--mapping_genomes: duplicate genome_id '${gid}' -- ids must be unique"
        }
        def fasta = file(f[1], checkIfExists: true)
        [ [ id: gid, taxon: (f.size() > 2 && f[2] != '-') ? f[2] : 'auto',
            large_index: fasta.size() > params.size_threshold_medium ], fasta ]
    }
    ch_genomes = Channel.fromList(genome_rows)

    INPUT_CHECK_WF(
        file(params.input, checkIfExists: true), params.raw_dir,
        params.strandedness, 'symbiont_mapping', params.sample
    )

    // Pre-aligned (BAM) rows are already committed to one genome and cannot
    // be scored against the others, so INPUT_CHECK_WF.out.reads silently
    // excludes them (see subworkflows/local/input_check). Silent is wrong
    // here -- say what got dropped rather than letting a sample vanish from
    // the comparison with no trace.
    INPUT_CHECK_WF.out.prealigned
        .map { meta, _bam -> meta.id }
        .collect()
        .subscribe { ids ->
            if( ids ) {
                log.warn "-entry symbiont_mapping: skipping ${ids.size()} pre-aligned " +
                          "sample(s) from --input (${ids.join(', ')}) -- a BAM is already " +
                          "aligned to one genome and cannot be compared across --mapping_genomes."
            }
        }

    MULTI_GENOME_MAPPING(INPUT_CHECK_WF.out.reads, ch_genomes)

    MULTI_GENOME_MAPPING.out.matrix.view      { "mapping matrix: ${it}" }
    MULTI_GENOME_MAPPING.out.best_genome.view { "best-genome calls: ${it}" }

    ch_versions = INPUT_CHECK_WF.out.versions
        .mix(MULTI_GENOME_MAPPING.out.versions)
        .unique()
        .collectFile(name: 'versions.yml', storeDir: "${params.tracedir}", sort: true)

    if( params.run_multiqc == null || params.run_multiqc ) {
        MULTIQC(
            MULTI_GENOME_MAPPING.out.multiqc.mix(ch_versions).collect().ifEmpty([]),
            file("${projectDir}/assets/multiqc_config.yml")
        )
    }
}

workflow.onComplete {
    def status = workflow.success ? 'completed' : 'FAILED'
    log.info """
    ${workflow.manifest.name} ${status}
      duration : ${workflow.duration}
      outdir   : ${params.outdir}
      strategy : ${params.outdir}/00_pipeline_info/strategy.yml
      versions : ${params.tracedir}/versions.yml
    """.stripIndent()
}
