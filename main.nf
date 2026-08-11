#!/usr/bin/env nextflow
/*
 * eukannot -- decision-driven eukaryote genome annotation + mRNA-seq pipeline
 * https://github.com/Ufungi/eukannot
 *
 * Entry points:
 *   (default)          full pipeline
 *   -entry strategy    resolve and print the execution plan, run nothing else
 *   -entry preflight   environment checks only
 */

nextflow.enable.dsl = 2

include { EUKANNOT     } from './workflows/eukannot'
include { GENOME_PREP  } from './subworkflows/local/genome_prep'
include { INPUT_CHECK_WF } from './subworkflows/local/input_check'

def helpMessage() {
    log.info """
    ${workflow.manifest.name} ${workflow.manifest.version}
    ${workflow.manifest.description}

    Usage:
      scripts/eukannot run . -profile singularity,local64 \\
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

    Profiles: singularity, docker, conda, gpu, local64, test, test_fungus
    """.stripIndent()
}

workflow {
    if( params.help ) {
        helpMessage()
        return
    }
    if( !(params.genome ?: params.fasta) ) {
        helpMessage()
        error "Missing required parameter: --genome (or --fasta)"
    }
    if( !params.input && !params.dry_run_strategy ) {
        helpMessage()
        error "Missing required parameter: --input (a samplesheet). " +
              "Use --dry_run_strategy to resolve a plan without one."
    }
    EUKANNOT()
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
