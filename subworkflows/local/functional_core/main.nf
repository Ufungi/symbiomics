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
include { SIGNALP6; TARGETP2; DEEPLOC2; DEEPTMHMM; PHOBIUS; EFFECTORP    } from '../../../modules/local/functional_predictors'
include { MEROPS_SCAN; CONSOLIDATE_PREDICTORS                            } from '../../../modules/local/functional_predictors'

def no(n) {
    return Channel.value(file("${projectDir}/assets/NO_FILE_${n}"))
}

workflow FUNCTIONAL_CORE {
    take:
    proteome   // path

    main:
    ch_versions = Channel.empty()
    wanted = params.functional_modules.tokenize(',')*.trim()

    ch_swissprot_hits = Channel.value(file("${projectDir}/assets/NO_FILE_SWISSPROT"))
    if( 'swissprot' in wanted ) {
        if( !params.swissprot_fasta || !file(params.swissprot_fasta).exists() ) {
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
        if( !params.eggnog_db_dir || !params.eggnog_dmnd || !file(params.eggnog_db_dir).exists() || !file(params.eggnog_dmnd).exists() ) {
            log.warn "[functional] eggnog requested but ${params.eggnog_db_dir} / ${params.eggnog_dmnd} not found -- SKIPPED"
        } else {
            EGGNOG_MAPPER(proteome, Channel.fromPath(params.eggnog_db_dir), Channel.fromPath(params.eggnog_dmnd))
            ch_versions = ch_versions.mix(EGGNOG_MAPPER.out.versions)
            ch_eggnog = EGGNOG_MAPPER.out.annotations
        }
    }

    ch_dbcan_overview = Channel.value(file("${projectDir}/assets/NO_FILE_DBCAN"))
    if( 'dbcan' in wanted ) {
        if( !params.cazy_db_dir || !file(params.cazy_db_dir).exists() ) {
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

    // Sequence-feature predictors (secretome / localisation / TM / effector /
    // peptidase). Licence-gated tools that are not installed yield an empty
    // file and are left out of the consensus; see docs/functional_annotation.md.
    ch_sp6 = no('SIGNALP6'); ch_tp2 = no('TARGETP2'); ch_dl2 = no('DEEPLOC2')
    ch_tmh = no('DEEPTMHMM'); ch_pho = no('PHOBIUS'); ch_eff = no('EFFECTORP'); ch_mer = no('MEROPS')
    if( 'signalp'   in wanted ) { SIGNALP6(proteome);   ch_sp6 = SIGNALP6.out.tsv;     ch_versions = ch_versions.mix(SIGNALP6.out.versions) }
    if( 'targetp'   in wanted ) { TARGETP2(proteome);   ch_tp2 = TARGETP2.out.tsv;     ch_versions = ch_versions.mix(TARGETP2.out.versions) }
    if( 'deeploc'   in wanted ) { DEEPLOC2(proteome);   ch_dl2 = DEEPLOC2.out.csv;     ch_versions = ch_versions.mix(DEEPLOC2.out.versions) }
    if( 'deeptmhmm' in wanted ) { DEEPTMHMM(proteome);  ch_tmh = DEEPTMHMM.out.topo;   ch_versions = ch_versions.mix(DEEPTMHMM.out.versions) }
    if( 'phobius'   in wanted ) { PHOBIUS(proteome);    ch_pho = PHOBIUS.out.tsv;      ch_versions = ch_versions.mix(PHOBIUS.out.versions) }
    if( 'effectorp' in wanted ) { EFFECTORP(proteome);  ch_eff = EFFECTORP.out.tsv;    ch_versions = ch_versions.mix(EFFECTORP.out.versions) }
    if( 'merops' in wanted ) {
        if( !params.merops_dmnd || !file(params.merops_dmnd).exists() ) {
            log.warn "[functional] merops requested but --merops_dmnd not found -- SKIPPED"
        } else {
            MEROPS_SCAN(proteome, Channel.fromPath(params.merops_dmnd))
            ch_mer = MEROPS_SCAN.out.tsv; ch_versions = ch_versions.mix(MEROPS_SCAN.out.versions)
        }
    }
    ch_localization = Channel.value(file("${projectDir}/assets/NO_FILE_LOCALIZATION"))
    if( wanted.any { it in ['signalp','targetp','deeploc','deeptmhmm','phobius','effectorp','merops'] } ) {
        CONSOLIDATE_PREDICTORS(proteome, ch_sp6, ch_tp2, ch_dl2, ch_tmh, ch_pho, ch_eff, ch_mer)
        ch_versions = ch_versions.mix(CONSOLIDATE_PREDICTORS.out.versions)
        ch_localization = CONSOLIDATE_PREDICTORS.out.tsv
    }

    // Phase B placeholders -- always empty until those modules exist.
    ch_funannotate2 = Channel.value(file("${projectDir}/assets/NO_FILE_FUNANNOTATE2"))
    ch_interpro     = Channel.value(file("${projectDir}/assets/NO_FILE_INTERPRO"))

    ASSIGN_PRODUCTS(proteome, ch_swissprot_hits, ch_eggnog, ch_funannotate2, ch_interpro, ch_dbcan_overview, ch_localization)
    ch_versions = ch_versions.mix(ASSIGN_PRODUCTS.out.versions)

    QC_PRODUCTS(ASSIGN_PRODUCTS.out.tsv)
    ch_versions = ch_versions.mix(QC_PRODUCTS.out.versions)

    RENDER_ANNOTATION_REPORT(ASSIGN_PRODUCTS.out.tsv)
    ch_versions = ch_versions.mix(RENDER_ANNOTATION_REPORT.out.versions)

    emit:
    localization    = ch_localization
    products        = ASSIGN_PRODUCTS.out.tsv
    products_qc     = QC_PRODUCTS.out.json
    annotation_report = RENDER_ANNOTATION_REPORT.out.html
    versions        = ch_versions
}
