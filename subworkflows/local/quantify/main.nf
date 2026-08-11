//
// Transcript assembly (StringTie) and counting (featureCounts / HTSeq).
//
// Both consume the single BAM channel from RNASEQ_ALIGN. Counting is the LAST
// consumer of that channel by design -- in the full pipeline it runs against
// the finished annotation, so it sits at the end of the DAG while the same
// BAMs have already served as annotation evidence upstream.
//

include { STRINGTIE_ASSEMBLE; STRINGTIE_MERGE } from '../../../modules/local/quantify'
include { INSPECT_ANNOTATION                  } from '../../../modules/local/quantify'
include { FEATURECOUNTS; HTSEQ_COUNT          } from '../../../modules/local/quantify'
include { MERGE_COUNTS                        } from '../../../modules/local/quantify'
include { MERGE_COUNTS as MERGE_COUNTS_HTSEQ  } from '../../../modules/local/quantify'

workflow RNASEQ_ASSEMBLE {
    take:
    bam         // [ meta, bam, index ]
    guide_gff   // path or NO_FILE

    main:
    ch_versions = Channel.empty()

    STRINGTIE_ASSEMBLE(bam, guide_gff)
    ch_versions = ch_versions.mix(STRINGTIE_ASSEMBLE.out.versions.first())

    STRINGTIE_MERGE(
        STRINGTIE_ASSEMBLE.out.gtf.map { _meta, gtf -> gtf }.collect(),
        guide_gff
    )
    ch_versions = ch_versions.mix(STRINGTIE_MERGE.out.versions)

    emit:
    per_sample = STRINGTIE_ASSEMBLE.out.gtf
    merged     = STRINGTIE_MERGE.out.gtf
    versions   = ch_versions
}

workflow QUANTIFY {
    take:
    bam             // [ meta, bam, index ]
    annotation      // path -- the annotation to count against
    strandedness    // collected *.strandedness.json (may be empty)

    main:
    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()
    ch_matrices = Channel.empty()

    // One row per sample, from the EFFECTIVE strandedness carried in meta --
    // which covers declared samples as well as probed ones. Folding the
    // inference JSONs in instead would leave declared samples blank.
    strand_tsv = bam
        .map { meta, _b, _i -> "${meta.id}\t${meta.strandedness}\n" }
        .collectFile(name: 'strandedness.tsv', sort: true)
        .ifEmpty(file("${projectDir}/assets/NO_FILE"))

    // Published GFF3s vary more than the counters assume. Rather than hardcode
    // -t exon -g gene_id and fail opaquely on a file that uses something else,
    // inspect the annotation and count on what it actually contains.
    INSPECT_ANNOTATION(annotation)
    ch_versions = ch_versions.mix(INSPECT_ANNOTATION.out.versions)

    spec = INSPECT_ANNOTATION.out.spec
        .map { json -> new groovy.json.JsonSlurperClassic().parseText(json.text) }
        .first()

    if( params.counter in ['featurecounts', 'both'] ) {
        FEATURECOUNTS(bam, annotation, spec)
        ch_versions = ch_versions.mix(FEATURECOUNTS.out.versions.first())
        ch_multiqc  = ch_multiqc.mix(FEATURECOUNTS.out.summary.map { _meta, s -> s })

        MERGE_COUNTS(
            FEATURECOUNTS.out.counts.map { _meta, counts -> counts }.collect(),
            'featurecounts',
            strand_tsv
        )
        ch_versions = ch_versions.mix(MERGE_COUNTS.out.versions)
        ch_matrices = ch_matrices.mix(MERGE_COUNTS.out.matrix)
    }

    if( params.counter in ['htseq', 'both'] ) {
        HTSEQ_COUNT(bam, annotation, spec)
        ch_versions = ch_versions.mix(HTSEQ_COUNT.out.versions.first())

        MERGE_COUNTS_HTSEQ(
            HTSEQ_COUNT.out.counts.map { _meta, counts -> counts }.collect(),
            'htseq',
            strand_tsv
        )
        ch_versions = ch_versions.mix(MERGE_COUNTS_HTSEQ.out.versions)
        ch_matrices = ch_matrices.mix(MERGE_COUNTS_HTSEQ.out.matrix)
    }

    emit:
    matrices = ch_matrices
    multiqc  = ch_multiqc
    versions = ch_versions
}
