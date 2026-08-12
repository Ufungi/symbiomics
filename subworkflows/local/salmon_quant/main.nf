//
// Salmon quantification -- an alternative to RNASEQ_ALIGN + QUANTIFY, not a
// modification of it. Deliberately independent: the whole point of the
// --quant_engine salmon path is fast quantification with no genome alignment
// and no BAM. It duplicates the (small) fastp/subsample logic that
// RNASEQ_ALIGN also uses rather than sharing state with it, so v0.1's tested
// HISAT2 path is never touched by this addition.
//
// transcript_fastas with length 1 -> ordinary single-haplotype quantification.
// transcript_fastas with length 2 -> diploid EM (docs/haplotypes.md option 2):
// IDs are checked disjoint, concatenated, and quantified as one index. The
// output here is the COMBINED gene matrix; splitting it into per-haplotype
// allele-specific counts is a separate step (bin/split_haplotype_counts.py)
// that needs the HA<->HB pairing table -- not yet available in this phase.
//

include { FASTQC; SEQKIT_SUBSAMPLE; FASTP } from '../../../modules/local/reads'
include { SALMON_CHECK_DISJOINT_IDS; SALMON_INDEX; SALMON_QUANT; SALMON_TO_GENE_MATRIX } from '../../../modules/local/salmon'
include { SPLIT_HAPLOTYPE_COUNTS as SPLIT_COUNTS } from '../../../modules/local/salmon'
include { SPLIT_HAPLOTYPE_COUNTS as SPLIT_TPM    } from '../../../modules/local/salmon'

workflow SALMON_QUANTIFY {
    take:
    reads              // [ meta, [fastq...] ]
    transcript_fastas  // list of 1 or 2 paths
    genome_fasta       // path, or NO_FILE if no decoy
    use_decoy          // val bool
    gff3               // path, or NO_FILE
    allele_pairs       // path, or NO_FILE -- from haplotype_pairing; only used when transcript_fastas.size() == 2

    main:
    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()

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

    // -- index (built once; a value channel so every sample's quant shares it)
    if( transcript_fastas.size() > 1 ) {
        SALMON_CHECK_DISJOINT_IDS(Channel.fromPath(transcript_fastas).collect())
        ch_versions = ch_versions.mix(SALMON_CHECK_DISJOINT_IDS.out.versions)
        ch_transcripts = SALMON_CHECK_DISJOINT_IDS.out.fasta
    } else {
        ch_transcripts = Channel.fromPath(transcript_fastas[0])
    }

    SALMON_INDEX(ch_transcripts, genome_fasta, use_decoy)
    ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)
    ch_index    = SALMON_INDEX.out.index.first()

    SALMON_QUANT(ch_trimmed, ch_index)
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())

    // -- gene-level matrix
    // Stage whole per-sample directories, not bare quant.sf files: every
    // sample's file is literally named "quant.sf", so staging the files
    // directly into one process work dir would collide on that basename.
    // A directory named after the sample id keeps them apart.
    quant_pairs = SALMON_QUANT.out.quant_dir
        .map { meta, dir -> "${meta.id}=${dir.name}/quant.sf" }
        .collect()
    quant_dirs = SALMON_QUANT.out.quant_dir.map { _meta, dir -> dir }.collect()

    SALMON_TO_GENE_MATRIX(quant_pairs, quant_dirs, gff3)
    ch_versions = ch_versions.mix(SALMON_TO_GENE_MATRIX.out.versions)

    // Diploid EM split (docs/haplotypes.md option 2): only meaningful with two
    // haplotypes' transcripts AND a pairing table -- otherwise the combined
    // matrix from above is already the final answer (a single haplotype has
    // nothing to split).
    ch_summed          = Channel.empty()
    ch_allele_specific = Channel.empty()
    ch_unpaired        = Channel.empty()
    if( transcript_fastas.size() > 1 && allele_pairs.name != 'NO_FILE' ) {
        SPLIT_COUNTS(SALMON_TO_GENE_MATRIX.out.counts, allele_pairs, 'counts')
        ch_versions = ch_versions.mix(SPLIT_COUNTS.out.versions)
        SPLIT_TPM(SALMON_TO_GENE_MATRIX.out.tpm, allele_pairs, 'tpm')
        ch_versions = ch_versions.mix(SPLIT_TPM.out.versions)

        ch_summed          = SPLIT_COUNTS.out.summed.mix(SPLIT_TPM.out.summed)
        ch_allele_specific = SPLIT_COUNTS.out.allele_specific.mix(SPLIT_TPM.out.allele_specific)
        ch_unpaired        = SPLIT_COUNTS.out.unpaired.mix(SPLIT_TPM.out.unpaired)
    }

    emit:
    counts          = SALMON_TO_GENE_MATRIX.out.counts        // combined HA+HB matrix
    tpm             = SALMON_TO_GENE_MATRIX.out.tpm
    summed          = ch_summed                                // PAV-inclusive total per allele pair
    allele_specific = ch_allele_specific                        // per-haplotype, labelled by pair
    unpaired        = ch_unpaired                                // genes with no partner (PAV candidates)
    quant_sf        = SALMON_QUANT.out.quant_sf
    multiqc         = ch_multiqc
    versions        = ch_versions
}
