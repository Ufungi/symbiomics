//
// The mRNA-seq arm: trim -> infer strandedness -> align -> sort -> index.
//
// The one structural property that matters: this subworkflow emits exactly ONE
// bam channel, which main.nf forks three ways (annotation evidence, StringTie,
// counting). Nothing downstream re-aligns. That is what "align once" means and
// it is why counting is deferred to the very end of the DAG rather than being
// done next to the alignment.
//
// Strandedness is a two-pass affair by necessity: --rna-strandness has to be
// set on the real alignment (it drives the XS tag StringTie reads), so pass 1
// aligns a small subsample unstranded purely to make the call.
//

include { FASTQC; SEQKIT_SUBSAMPLE; FASTP        } from '../../../modules/local/reads'
include { HISAT2_BUILD; HISAT2_ALIGN             } from '../../../modules/local/align'
include { SAMTOOLS_SORT; SAMTOOLS_INDEX          } from '../../../modules/local/align'
include { BAM_STATS; INFER_STRANDEDNESS          } from '../../../modules/local/align'
include { HISAT2_ALIGN as HISAT2_ALIGN_PROBE     } from '../../../modules/local/align'
include { SEQKIT_SUBSAMPLE as SEQKIT_PROBE       } from '../../../modules/local/reads'

workflow RNASEQ_ALIGN {
    take:
    reads          // [ meta, [fastq...] ]
    genome         // [ fasta, fai ]
    strategy       // value channel carrying the decoded strategy.json Map
    prealigned     // [ meta, bam ]

    main:
    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()

    // The plan is produced by a process, so it can only be read through a
    // channel. Derive each decision as its own value channel and feed it in as
    // a process `val` input -- that is what makes the policy engine actually
    // govern execution rather than just describe it.
    use_mmap    = strategy.map { plan -> plan?.decisions?.aligner?.mmap_index ?: false }
    large_index = strategy.map { plan -> plan?.decisions?.aligner?.large_index ?: false }
    bam_index   = strategy.map { plan -> plan?.decisions?.bam_index?.choice ?: 'csi' }

    // ------------------------------------------------------------- 30 reads
    FASTQC(reads)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    ch_multiqc  = ch_multiqc.mix(FASTQC.out.zip.map { _meta, zip -> zip })

    ch_input = reads
    if( params.test_mode ) {
        SEQKIT_SUBSAMPLE(reads, params.test_reads, params.test_sampler, params.test_seed)
        ch_versions = ch_versions.mix(SEQKIT_SUBSAMPLE.out.versions.first())
        ch_input = SEQKIT_SUBSAMPLE.out.reads
    }

    if( params.trimmer == 'fastp' ) {
        FASTP(ch_input)
        ch_trimmed  = FASTP.out.reads
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
        ch_multiqc  = ch_multiqc.mix(FASTP.out.json.map { _meta, json -> json })
    } else {
        ch_trimmed = ch_input
    }

    // -------------------------------------------------------------- 40 index
    if( params.hisat2_index ) {
        // Reusing a prebuilt index is the difference between minutes and days
        // on a multi-gigabase genome. --hisat2_index is a PREFIX in hisat2's own
        // sense (the thing before `.1.ht2`), so split it into the directory to
        // stage and the basename to pass to -x.
        def idx_path   = file(params.hisat2_index)
        def idx_dir    = idx_path.isDirectory() ? idx_path : idx_path.parent
        def idx_prefix = idx_path.isDirectory() ? 'genome' : idx_path.name
        def found = idx_dir.list().findAll { it.startsWith(idx_prefix) && (it.endsWith('.ht2') || it.endsWith('.ht2l')) }
        if( !found ) {
            error("""
            --hisat2_index ${params.hisat2_index} does not resolve to a HISAT2 index.
              looked in : ${idx_dir}
              for       : ${idx_prefix}.*.ht2 / .ht2l
            Give the prefix, not a single .ht2 file (e.g. /path/genome, not /path/genome.1.ht2).
            """.stripIndent())
        }
        log.info "[align] reusing HISAT2 index: ${idx_dir}/${idx_prefix} (${found.size()} files)"
        ch_index  = Channel.value(idx_dir)
        ch_prefix = Channel.value(idx_prefix)
    } else {
        HISAT2_BUILD(genome.map { fasta, _fai -> fasta }, large_index)
        ch_index    = HISAT2_BUILD.out.index
        ch_prefix   = Channel.value('genome')
        ch_versions = ch_versions.mix(HISAT2_BUILD.out.versions)
    }

    // ------------------------------------------------- 35 strandedness probe
    needs_probe = ch_trimmed.filter { meta, _fq -> meta.strandedness == 'auto' }

    SEQKIT_PROBE(needs_probe, params.strandedness_subsample, 'head', params.test_seed)
    ch_versions = ch_versions.mix(SEQKIT_PROBE.out.versions.first())

    // Pass 1 is deliberately unstranded: the whole point is to measure, not
    // assume. The declared value is preserved under a separate key -- if the
    // override leaked into `strandedness` the probe would compare its own
    // measurement against itself and no declared/inferred conflict could ever
    // be detected.
    probe_reads = SEQKIT_PROBE.out.reads
        .map { meta, fq ->
            [ meta + [strandedness: 'unstranded', declared: meta.strandedness], fq ]
        }

    HISAT2_ALIGN_PROBE(probe_reads, ch_index, ch_prefix, use_mmap)
    ch_versions = ch_versions.mix(HISAT2_ALIGN_PROBE.out.versions.first())

    INFER_STRANDEDNESS(HISAT2_ALIGN_PROBE.out.bam, genome)
    ch_versions = ch_versions.mix(INFER_STRANDEDNESS.out.versions.first())
    ch_multiqc  = ch_multiqc.mix(INFER_STRANDEDNESS.out.json.map { _meta, json -> json })

    resolved_strand = INFER_STRANDEDNESS.out.json
        .map { meta, json ->
            def call = new groovy.json.JsonSlurper().parse(json)
            [ meta.id, call.effective ]
        }

    // Rows that declared a strandedness keep it; probed rows take the call.
    declared = ch_trimmed
        .filter { meta, _fq -> meta.strandedness != 'auto' }
        .map { meta, fq -> [ meta, fq ] }

    probed = ch_trimmed
        .filter { meta, _fq -> meta.strandedness == 'auto' }
        .map { meta, fq -> [ meta.id, meta, fq ] }
        .join(resolved_strand, by: 0)
        .map { _id, meta, fq, call -> [ meta + [strandedness: call], fq ] }

    ch_for_align = declared.mix(probed)

    // ------------------------------------------------------------- 40 align
    HISAT2_ALIGN(ch_for_align, ch_index, ch_prefix, use_mmap)
    ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions.first())
    ch_multiqc  = ch_multiqc.mix(HISAT2_ALIGN.out.summary.map { _meta, txt -> txt })

    SAMTOOLS_SORT(HISAT2_ALIGN.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    // Pre-aligned rows rejoin here, having skipped trim and align entirely.
    ch_all_bam = SAMTOOLS_SORT.out.bam.mix(prealigned)

    SAMTOOLS_INDEX(ch_all_bam, bam_index)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    BAM_STATS(SAMTOOLS_INDEX.out.bam)
    ch_versions = ch_versions.mix(BAM_STATS.out.versions.first())
    ch_multiqc  = ch_multiqc.mix(BAM_STATS.out.stats, BAM_STATS.out.flagstat, BAM_STATS.out.idxstats)

    emit:
    bam          = SAMTOOLS_INDEX.out.bam         // [ meta, bam, index ]  <- THE channel
    strandedness = INFER_STRANDEDNESS.out.json.map { _meta, json -> json }
    multiqc      = ch_multiqc
    versions     = ch_versions
}
