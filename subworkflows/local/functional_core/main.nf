//
// Shared functional-annotation core -- the taxon != fungi branch of
// -entry functional, and also the modules a fungal run could layer on top of
// (Phase A only for now; see docs/functional_annotation.md for the Phase B
// modules -- InterProScan6, KofamScan, dbCAN v5 -- gated behind the same
// `functional_modules` list once their databases are provisioned).
//
// Every module here is independently optional (`params.functional_modules`),
// and a module whose database is not provisioned is SKIPPED, not a hard
// failure, exactly like the on_missing_tool convention from the mRNA-seq arm.
//

include { DIAMOND_MAKEDB; DIAMOND_BLASTP_SWISSPROT; PARSE_SWISSPROT_HITS } from '../../../modules/local/functional'
include { EGGNOG_MAPPER; RUN_DBCAN; ASSIGN_PRODUCTS; QC_PRODUCTS         } from '../../../modules/local/functional'
include { RENDER_ANNOTATION_REPORT                                      } from '../../../modules/local/functional'

workflow FUNCTIONAL_CORE {
    take:
    proteome   // path

    main:
    ch_versions = Channel.empty()
    wanted = params.functional_modules.tokenize(',')*.trim()

    ch_swissprot_hits = Channel.value(file("${projectDir}/assets/NO_FILE_SWISSPROT"))
    if( 'swissprot' in wanted ) {
        if( !file(params.swissprot_fasta).exists() ) {
            log.warn "[functional] swissprot requested but ${params.swissprot_fasta} does not exist -- SKIPPED"
        } else {
            DIAMOND_MAKEDB(Channel.fromPath(params.swissprot_fasta))
            ch_versions = ch_versions.mix(DIAMOND_MAKEDB.out.versions)

            DIAMOND_BLASTP_SWISSPROT(proteome, DIAMOND_MAKEDB.out.db)
            ch_versions = ch_versions.mix(DIAMOND_BLASTP_SWISSPROT.out.versions)

            PARSE_SWISSPROT_HITS(DIAMOND_BLASTP_SWISSPROT.out.tsv, Channel.fromPath(params.swissprot_fasta))
            ch_versions = ch_versions.mix(PARSE_SWISSPROT_HITS.out.versions)
            ch_swissprot_hits = PARSE_SWISSPROT_HITS.out.tsv
        }
    }

    ch_eggnog = Channel.value(file("${projectDir}/assets/NO_FILE_EGGNOG"))
    if( 'eggnog' in wanted ) {
        if( !file(params.eggnog_db_dir).exists() || !file(params.eggnog_dmnd).exists() ) {
            log.warn "[functional] eggnog requested but ${params.eggnog_db_dir} / ${params.eggnog_dmnd} not found -- SKIPPED"
        } else {
            EGGNOG_MAPPER(proteome, Channel.fromPath(params.eggnog_db_dir), Channel.fromPath(params.eggnog_dmnd))
            ch_versions = ch_versions.mix(EGGNOG_MAPPER.out.versions)
            ch_eggnog = EGGNOG_MAPPER.out.annotations
        }
    }

    ch_dbcan_overview = Channel.value(file("${projectDir}/assets/NO_FILE_DBCAN"))
    if( 'dbcan' in wanted ) {
        if( !file(params.cazy_db_dir).exists() ) {
            log.warn "[functional] dbcan requested but ${params.cazy_db_dir} not found -- SKIPPED"
        } else {
            RUN_DBCAN(proteome, Channel.fromPath(params.cazy_db_dir))
            ch_versions = ch_versions.mix(RUN_DBCAN.out.versions)
            // overview.txt is emitted `optional: true` -- fall back to the
            // sentinel if a run produces no CAZyme calls at all, same as
            // every other optional functional module here.
            ch_dbcan_overview = RUN_DBCAN.out.overview.ifEmpty(file("${projectDir}/assets/NO_FILE_DBCAN"))
        }
    }

    // Phase B placeholders -- always empty until those modules exist.
    ch_funannotate2 = Channel.value(file("${projectDir}/assets/NO_FILE_FUNANNOTATE2"))
    ch_interpro     = Channel.value(file("${projectDir}/assets/NO_FILE_INTERPRO"))

    ASSIGN_PRODUCTS(proteome, ch_swissprot_hits, ch_eggnog, ch_funannotate2, ch_interpro, ch_dbcan_overview)
    ch_versions = ch_versions.mix(ASSIGN_PRODUCTS.out.versions)

    QC_PRODUCTS(ASSIGN_PRODUCTS.out.tsv)
    ch_versions = ch_versions.mix(QC_PRODUCTS.out.versions)

    RENDER_ANNOTATION_REPORT(ASSIGN_PRODUCTS.out.tsv)
    ch_versions = ch_versions.mix(RENDER_ANNOTATION_REPORT.out.versions)

    emit:
    products        = ASSIGN_PRODUCTS.out.tsv
    products_qc     = QC_PRODUCTS.out.json
    annotation_report = RENDER_ANNOTATION_REPORT.out.html
    versions        = ch_versions
}
