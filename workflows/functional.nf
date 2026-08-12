//
// -entry functional -- protein-in functional annotation, decoupled entirely
// from structural annotation. No genome, no samplesheet, no BRAKER: if a good
// gene model already exists (as it does for the reference Pinus densiflora
// assembly this pipeline was built against), there is no reason to run
// structural annotation just to get back to the protein set you already have.
//
// taxon=fungi routes to funannotate2 + funannotate2-addons (F4/Phase B --
// wraps eggNOG-mapper and InterProScan itself, so FUNCTIONAL_CORE must not
// duplicate them there). Not yet implemented; errors clearly rather than
// silently running the wrong branch.
//
// taxon != fungi routes to FUNCTIONAL_CORE: DIAMOND-vs-Swiss-Prot,
// eggNOG-mapper, dbCAN today (Phase A, zero new downloads); InterProScan6 /
// KofamScan / dbCAN v5 once their databases are provisioned (Phase B).
//

include { FUNCTIONAL_CORE } from '../subworkflows/local/functional_core'
include { MULTIQC          } from '../modules/local/quantify'

workflow FUNCTIONAL {
    if( !params.proteome ) {
        error "Missing required parameter: --proteome <protein.fasta>"
    }
    proteome_file = file(params.proteome, checkIfExists: true)

    def taxon = params.taxon ?: 'other'
    if( taxon == 'fungi' ) {
        error """
        --taxon fungi routes to funannotate2 + funannotate2-addons, which is not
        yet implemented (planned: F4 / Phase B, see docs/roadmap.md).
        For now, run with --taxon plant/animal/other to use the shared
        DIAMOND-vs-Swiss-Prot + eggNOG-mapper + dbCAN core instead.
        """.stripIndent()
    }

    FUNCTIONAL_CORE(Channel.fromPath(proteome_file))

    ch_versions = FUNCTIONAL_CORE.out.versions
        .unique()
        .collectFile(name: 'versions.yml', storeDir: "${params.tracedir}", sort: true)

    if( params.run_multiqc == null || params.run_multiqc ) {
        MULTIQC(
            FUNCTIONAL_CORE.out.products
                .mix(FUNCTIONAL_CORE.out.products_qc)
                .mix(ch_versions)
                .collect().ifEmpty([]),
            file("${projectDir}/assets/multiqc_config.yml")
        )
    }

    FUNCTIONAL_CORE.out.products_qc.view { qc ->
        "functional annotation done -- QC summary: ${qc}"
    }
}
