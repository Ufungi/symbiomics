//
// Map one set of transcriptome reads against N candidate genomes.
//
// This is the "symbiomics" part of symbiomics: a matsutake/pine-style
// dual-organism sample (or a set of candidate symbiont genomes of unknown
// identity) needs the SAME reads scored against every candidate, not one
// genome chosen up front. trim -> build N indices -> align every sample
// against every genome -> per-sample/per-genome mapping-rate matrix and a
// best-genome call.
//
// Deliberately NOT folded into RNASEQ_ALIGN: that subworkflow's whole
// contract is "align once" against a single genome, forked downstream to
// three consumers. Here the fork happens on the INPUT side instead (one
// sample x many genomes), which is a different shape, not a variant of the
// same one.
//

include { FASTP                            } from '../../../modules/local/reads'
include { HISAT2_BUILD_MULTI               } from '../../../modules/local/multi_genome_map'
include { HISAT2_ALIGN_MULTI               } from '../../../modules/local/multi_genome_map'
include { SUMMARIZE_MULTI_GENOME_MAPPING   } from '../../../modules/local/multi_genome_map'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_MULTI   } from '../../../modules/local/align'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MULTI } from '../../../modules/local/align'
include { BAM_STATS as BAM_STATS_MULTI           } from '../../../modules/local/align'

workflow MULTI_GENOME_MAPPING {
    take:
    reads    // [ meta, [fastq...] ]                     -- INPUT_CHECK_WF.out.reads
    genomes  // [ genome_meta, fasta ]                    -- genome_meta: [id, taxon, large_index]

    main:
    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()

    if( params.trimmer == 'fastp' ) {
        FASTP(reads)
        ch_trimmed  = FASTP.out.reads
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
        ch_multiqc  = ch_multiqc.mix(FASTP.out.json.map { _meta, json -> json })
    } else {
        ch_trimmed = reads
    }

    // ------------------------------------------------------------- N indices
    HISAT2_BUILD_MULTI(genomes)
    ch_versions = ch_versions.mix(HISAT2_BUILD_MULTI.out.versions.first())

    // ------------------------------------------------- every sample x genome
    // .combine() cross-joins: with G genomes this fans one sample out into G
    // alignment tasks, which is the entire point -- the same reads scored
    // against every candidate rather than committed to one up front.
    combo = ch_trimmed
        .combine(HISAT2_BUILD_MULTI.out.index)
        .map { meta, fq, genome_meta, index ->
            def combo_meta = meta + [
                id         : "${meta.id}__${genome_meta.id}",
                sample_id  : meta.id,
                genome_id  : genome_meta.id,
                large_index: genome_meta.large_index,
            ]
            [ combo_meta, fq, index ]
        }

    HISAT2_ALIGN_MULTI(combo)
    ch_versions = ch_versions.mix(HISAT2_ALIGN_MULTI.out.versions.first())
    ch_multiqc  = ch_multiqc.mix(HISAT2_ALIGN_MULTI.out.summary.map { _meta, txt -> txt })

    SAMTOOLS_SORT_MULTI(HISAT2_ALIGN_MULTI.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_MULTI.out.versions.first())

    bam_index = params.bam_index in ['csi', 'bai', 'both'] ? params.bam_index : 'csi'
    SAMTOOLS_INDEX_MULTI(SAMTOOLS_SORT_MULTI.out.bam, Channel.value(bam_index))
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_MULTI.out.versions.first())

    BAM_STATS_MULTI(SAMTOOLS_INDEX_MULTI.out.bam)
    ch_versions = ch_versions.mix(BAM_STATS_MULTI.out.versions.first())
    ch_multiqc  = ch_multiqc.mix(BAM_STATS_MULTI.out.stats, BAM_STATS_MULTI.out.flagstat, BAM_STATS_MULTI.out.idxstats)

    // ----------------------------------------------------------- comparison
    manifest = HISAT2_ALIGN_MULTI.out.summary
        .map { meta, summary -> "${meta.sample_id}\t${meta.genome_id}\t${summary.name}\n" }
        .collectFile(name: 'manifest.tsv', sort: true, newLine: false,
                     seed: "sample_id\tgenome_id\tsummary_file\n")

    SUMMARIZE_MULTI_GENOME_MAPPING(
        HISAT2_ALIGN_MULTI.out.summary.map { _meta, summary -> summary }.collect(),
        manifest
    )
    ch_versions = ch_versions.mix(SUMMARIZE_MULTI_GENOME_MAPPING.out.versions)

    emit:
    bam          = SAMTOOLS_INDEX_MULTI.out.bam     // [ combo_meta, bam, index ]
    matrix       = SUMMARIZE_MULTI_GENOME_MAPPING.out.matrix
    best_genome  = SUMMARIZE_MULTI_GENOME_MAPPING.out.best_genome
    summary_json = SUMMARIZE_MULTI_GENOME_MAPPING.out.json
    multiqc      = ch_multiqc
    versions     = ch_versions
}
