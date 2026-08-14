//
// Functional annotation, Phase A: DIAMOND-vs-Swiss-Prot, eggNOG-mapper,
// dbCAN -- every one of these runs against a database already on this lab's
// disk, so `-entry functional` is usable the day this ships, before any of
// Phase B's ~80 GB of new downloads (InterProScan6, dbCAN v5, KofamScan,
// funannotate2) land. See docs/functional_annotation.md.
//

process DIAMOND_MAKEDB {
    label 'process_medium'
    tag "${fasta.name}"

    input:
    path fasta

    output:
    path '*.dmnd'      , emit: db
    path 'versions.yml', emit: versions

    script:
    def base = fasta.baseName
    """
    diamond makedb --in '${fasta}' -d ${base} --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond version | sed 's/diamond version //')
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta.baseName}.dmnd
    echo '"${task.process}": {diamond: stub}' > versions.yml
    """
}

process DIAMOND_BLASTP_SWISSPROT {
    label 'process_medium'
    tag "${proteome.name}"

    input:
    path proteome
    path swissprot_dmnd

    output:
    path 'swissprot.diamond.tsv', emit: tsv
    path 'versions.yml'         , emit: versions

    script:
    """
    diamond blastp \\
        -q '${proteome}' -d '${swissprot_dmnd}' \\
        -o swissprot.diamond.tsv \\
        --outfmt 6 qseqid sseqid pident length evalue bitscore \\
        --max-target-seqs 1 --evalue 1e-5 --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond version | sed 's/diamond version //')
    END_VERSIONS
    """

    stub:
    """
    touch swissprot.diamond.tsv
    echo '"${task.process}": {diamond: stub}' > versions.yml
    """
}

process PARSE_SWISSPROT_HITS {
    label 'process_single'
    tag "${diamond_tsv.name}"

    input:
    path diamond_tsv
    path swissprot_fasta

    output:
    path 'swissprot_hits.tsv', emit: tsv
    path 'versions.yml'      , emit: versions

    script:
    """
    parse_diamond_swissprot.py \\
        --diamond-tsv '${diamond_tsv}' --swissprot-fasta '${swissprot_fasta}' \\
        --out swissprot_hits.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf 'query_id\\tswissprot_accession\\tswissprot_entry\\tpident\\tevalue\\tbitscore\\tdescription\\torganism\\n' > swissprot_hits.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}

process EGGNOG_MAPPER {
    label 'process_high'
    tag "${proteome.name}"

    input:
    path proteome
    path eggnog_db_dir
    path eggnog_dmnd

    output:
    path '*.emapper.annotations', emit: annotations
    path '*.emapper.seed_orthologs', emit: seed_orthologs, optional: true
    path 'versions.yml', emit: versions

    script:
    """
    emapper.py -m diamond --itype proteins \\
        -i '${proteome}' -o symbiomics --cpu ${task.cpus} \\
        --dmnd_db '${eggnog_dmnd}' --data_dir '${eggnog_db_dir}' \\
        --evalue 1e-5 --temp_dir \$PWD

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: \$(emapper.py --version 2>&1 | head -1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
    END_VERSIONS
    """

    stub:
    """
    printf '## stub\\n#query\\tseed_ortholog\\tevalue\\tscore\\teggNOG_OGs\\tmax_annot_lvl\\tCOG_category\\tDescription\\tPreferred_name\\tGOs\\tEC\\tKEGG_ko\\n' > symbiomics.emapper.annotations
    echo '"${task.process}": {eggnog-mapper: stub}' > versions.yml
    """
}

process RUN_DBCAN {
    label 'process_medium'
    tag "${proteome.name}"
    // Phase A uses the old (v3-line) run_dbcan CLI against this lab's existing
    // 2021-vintage CAZy database. dbCAN's CLI and DB layout changed materially
    // in v4/v5 (see docs/functional_annotation.md) -- do not point this
    // container at a newer DB or vice versa.
    //
    // All three callers requested (hmmer + diamond + eCAMI), not hmmer alone --
    // this is what dbCAN3 is actually designed around: overview.txt's whole
    // point is a #ofTools consensus column across all three, not a single
    // caller's raw output. Verified directly against the pinned container
    // (dbcan:3.0.7--pyh5e36f6f_0): eCAMI's default kmer database ("CAZyme")
    // ships inside the eCAMI package itself, not something this lab has to
    // provision separately -- a real run against this DB directory completed
    // in ~35s with all three columns populated, no missing-input error.

    input:
    path proteome
    path cazy_db_dir

    output:
    path 'dbcan_out/overview.txt', emit: overview, optional: true
    path 'dbcan_out'             , emit: outdir
    path 'versions.yml'          , emit: versions

    script:
    """
    run_dbcan '${proteome}' protein \\
        --db_dir '${cazy_db_dir}' --out_dir dbcan_out \\
        --tools hmmer diamond eCAMI \\
        --hmm_cpu ${task.cpus} --dia_cpu ${task.cpus} --eCAMI_jobs ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        run_dbcan: 3.0.7 (Phase A, hmmer+diamond+eCAMI consensus, pinned to this lab's 2021-vintage CAZy DB)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p dbcan_out
    printf 'Gene ID\\tHMMER\\n' > dbcan_out/overview.txt
    echo '"${task.process}": {run_dbcan: stub}' > versions.yml
    """
}

process ASSIGN_PRODUCTS {
    label 'process_single'
    tag "${proteome.name}"

    input:
    path proteome
    path swissprot_hits   // NO_FILE if not run
    path eggnog_annotations // NO_FILE if not run
    path funannotate2_tsv  // NO_FILE if not run (fungal branch, F4)
    path interpro_tsv      // NO_FILE if not run (Phase B)
    path dbcan_overview    // NO_FILE if not run

    output:
    path 'products.tsv', emit: tsv
    path 'versions.yml', emit: versions

    script:
    // Five independently-optional inputs can all be "not provided" at once,
    // and staging multiple files literally named the same thing into one task
    // directory collides -- so each has its OWN distinctly-named sentinel
    // (assets/NO_FILE_*) rather than sharing the single NO_FILE used
    // elsewhere in this pipeline. Presence is tested by size, not name.
    def sp = swissprot_hits.size() > 0 ? "--swissprot '${swissprot_hits}'" : ''
    def eg = eggnog_annotations.size() > 0 ? "--eggnog '${eggnog_annotations}'" : ''
    def f2 = funannotate2_tsv.size() > 0 ? "--funannotate2 '${funannotate2_tsv}'" : ''
    def ip = interpro_tsv.size() > 0 ? "--interpro '${interpro_tsv}'" : ''
    def dc = dbcan_overview.size() > 0 ? "--dbcan '${dbcan_overview}'" : ''
    """
    assign_products.py --proteome '${proteome}' ${sp} ${eg} ${f2} ${ip} ${dc} --out products.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf 'protein_id\\tproduct\\tsource\\tgo_terms\\tec_number\\tkegg_ko\\tswissprot_hit\\tswissprot_pident\\tdbcan_family\\tdbcan_tools\\tsources_with_evidence\\tn_sources\\n' > products.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}

process RENDER_ANNOTATION_REPORT {
    label 'process_single'
    tag "${products.name}"

    input:
    path products

    output:
    path 'annotation_report.html', emit: html
    path 'versions.yml'          , emit: versions

    script:
    """
    render_annotation_report.py --products '${products}' --out annotation_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch annotation_report.html
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}

process QC_PRODUCTS {
    label 'process_single'
    tag "${products.name}"

    input:
    path products

    output:
    path 'products_qc.json', emit: json
    path 'versions.yml'     , emit: versions

    script:
    """
    qc_products.py --products '${products}' \\
        --max-hypothetical-frac ${params.max_hypothetical_frac} \\
        --out products_qc.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo '{"total_proteins":0,"by_source":{},"hypothetical_fraction":0.0,"over_threshold":false}' > products_qc.json
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
