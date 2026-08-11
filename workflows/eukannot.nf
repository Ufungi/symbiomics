//
// The main workflow.
//
// Shape of the v0.1 DAG:
//
//   samplesheet -> INPUT_CHECK ---------------------------+
//                                                          |
//   genome ------> GENOME_PREP -> DECIDE_STRATEGY -> strategy.yml
//                                                          |
//                        (strategy drives every choice below)
//                                                          |
//   reads -> trim -> strandedness probe -> HISAT2 -> sort -> INDEX
//                                                          |
//                       +---------------- BAM -------------+
//                       |                 |                |
//                consumer 1          consumer 2       consumer 3
//             (annotation evidence,  StringTie      counting, held
//              M3 -- not in v0.1)                   to the end)
//
// Structural annotation (M3+) plugs in between StringTie and QUANTIFY without
// changing the BAM channel, which is why the fork exists already.
//

include { INPUT_CHECK_WF          } from '../subworkflows/local/input_check'
include { GENOME_PREP             } from '../subworkflows/local/genome_prep'
include { RNASEQ_ALIGN            } from '../subworkflows/local/rnaseq_align'
include { RNASEQ_ASSEMBLE         } from '../subworkflows/local/quantify'
include { QUANTIFY                } from '../subworkflows/local/quantify'
include { MULTIQC                 } from '../modules/local/quantify'

def enabled(String step) {
    def wanted  = params.steps == 'all' ? null : params.steps.tokenize(',')*.trim()
    def skipped = params.skip_steps ? params.skip_steps.tokenize(',')*.trim() : []
    if( step in skipped ) return false
    return wanted == null || step in wanted
}

workflow EUKANNOT {

    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()

    // Hardlink publishing silently fails across filesystems; catch it here
    // rather than after the first 40 GB BAM has been written.
    if( params.publish_mode == 'link' ) {
        def work_dev = java.nio.file.Files.getFileStore(workflow.workDir).name()
        def out_dir  = file(params.outdir); out_dir.mkdirs()
        def out_dev  = java.nio.file.Files.getFileStore(out_dir.toPath()).name()
        if( work_dev != out_dev ) {
            error """
            --publish_mode link needs the work directory and --outdir on one filesystem.
              work   : ${workflow.workDir}  (${work_dev})
              outdir : ${out_dir}  (${out_dev})
            fix: export NXF_WORK=${out_dir}/work   -- or use --publish_mode copy
            """.stripIndent()
        }
    }

    // ------------------------------------------------------------- 00 input
    genome_file = file(params.genome ?: params.fasta, checkIfExists: true)
    genome_id   = params.genome_id ?: genome_file.simpleName
    genome_meta = [ id: genome_id, species: params.species ?: genome_id ]

    INPUT_CHECK_WF(
        file(params.input, checkIfExists: true),
        params.raw_dir,
        params.strandedness,
        genome_id,
        params.sample
    )
    ch_versions = ch_versions.mix(INPUT_CHECK_WF.out.versions)

    // Evidence available to the policy engine. `n_samples` is resolved lazily
    // from the validated sheet rather than guessed from the raw file.
    evidence = INPUT_CHECK_WF.out.resolved
        .map { tsv ->
            def rows = tsv.readLines().findAll { it && !it.startsWith('sample_id') }
            [
                has_rna    : rows.size() > 0,
                has_protein: (params.proteins ?: '') != '' && params.proteins != 'auto',
                has_isoseq : rows.any { it.split('\t').size() > 7 && it.split('\t')[7] != 'short' },
                n_samples  : rows.size(),
            ]
        }

    // Byte size of a reused index feeds the --mm decision.
    index_bytes = params.hisat2_index
        ? file(params.hisat2_index).parent.listFiles()
              .findAll { it.name.endsWith('.ht2') || it.name.endsWith('.ht2l') }
              .sum { it.size() } ?: 0
        : 0

    // ------------------------------------------------ 10 genome + strategy
    GENOME_PREP(
        Channel.value([ genome_meta, genome_file ]),
        evidence,
        index_bytes
    )
    ch_versions = ch_versions.mix(GENOME_PREP.out.versions)

    // A value channel: the genome is read by several processes and a queue
    // channel would be drained by the first of them.
    ch_genome = GENOME_PREP.out.fai.map { fai -> [ genome_file, fai ] }.first()

    if( params.dry_run_strategy ) {
        GENOME_PREP.out.yml.view { "strategy written: ${it}" }
        return
    }

    // The decoded strategy is a VALUE channel so every consumer sees the same
    // map without re-deriving decisions and without consuming it.
    ch_strategy = GENOME_PREP.out.strategy.map { _meta, plan -> plan }.first()

    // ---------------------------------------------------- 30-40 mRNA-seq arm
    if( enabled('align') ) {
        RNASEQ_ALIGN(
            INPUT_CHECK_WF.out.reads,
            ch_genome,
            ch_strategy,
            INPUT_CHECK_WF.out.prealigned
        )
        ch_versions = ch_versions.mix(RNASEQ_ALIGN.out.versions)
        ch_multiqc  = ch_multiqc.mix(RNASEQ_ALIGN.out.multiqc)

        // ------------------------------------------------- THE three-way fork
        RNASEQ_ALIGN.out.bam
            .multiMap { meta, bam, index ->
                to_evidence: [ meta, bam, index ]   // consumer 1 (M3: BRAKER)
                to_assemble: [ meta, bam, index ]   // consumer 2
                to_quant   : [ meta, bam, index ]   // consumer 3, held to the end
            }
            .set { bam_fork }

        guide = params.reference_gff
            ? file(params.reference_gff, checkIfExists: true)
            : file("${projectDir}/assets/NO_FILE")

        // ------------------------------------------------------- 45 assemble
        if( enabled('assemble') && params.run_stringtie ) {
            RNASEQ_ASSEMBLE(
                bam_fork.to_assemble.filter { meta, _b, _i -> 'stringtie' in meta.use_for },
                guide
            )
            ch_versions = ch_versions.mix(RNASEQ_ASSEMBLE.out.versions)
            ch_annotation = RNASEQ_ASSEMBLE.out.merged
        } else {
            ch_annotation = Channel.empty()
        }

        // ------------------------------------------------------- 95 quantify
        // In v0.1 the countable annotation is either the user's reference GFF
        // or the StringTie merge. From M3 it becomes the finished gene set.
        quant_gff = params.quant_gff        ? Channel.value(file(params.quant_gff, checkIfExists: true))
                  : params.reference_gff    ? Channel.value(file(params.reference_gff, checkIfExists: true))
                  : ch_annotation

        if( enabled('quantify') && params.counter != 'none' ) {
            QUANTIFY(
                bam_fork.to_quant.filter { meta, _b, _i -> 'count' in meta.use_for },
                quant_gff,
                RNASEQ_ALIGN.out.strandedness
            )
            ch_versions = ch_versions.mix(QUANTIFY.out.versions)
            ch_multiqc  = ch_multiqc.mix(QUANTIFY.out.multiqc)
        }
    }

    // --------------------------------------------------------------- 99 report
    ch_versions
        .unique()
        .collectFile(name: 'versions.yml', storeDir: "${params.tracedir}", sort: true)
        .set { ch_versions_file }

    if( params.run_multiqc == null || params.run_multiqc ) {
        MULTIQC(
            ch_multiqc.mix(ch_versions_file).collect().ifEmpty([]),
            file("${projectDir}/assets/multiqc_config.yml")
        )
    }
}
