//
// -entry download_dbs -- provision functional-annotation databases.
// Select with --download_dbs swissprot,eggnog (default: swissprot, ~90 MB).
//   swissprot  UniProt Swiss-Prot FASTA (curl)
//   eggnog     eggNOG-mapper DB via its own download_eggnog_data.py (tens of GB!)
// dbCAN/CAZy, SignalP/TargetP/DeepLoc/DeepTMHMM/Phobius are licence- or
// version-specific and are NOT downloaded here; see docs/functional_annotation.md.
//

process DOWNLOAD_SWISSPROT {
    tag 'swissprot'
    storeDir "${params.db_dir}/swissprot"
    output: path 'uniprot_sprot.fasta'
    script:
    """
    curl -fL --retry 5 -o uniprot_sprot.fasta.gz https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
    gunzip uniprot_sprot.fasta.gz
    """
    stub:
    "echo '>sp|stub' > uniprot_sprot.fasta"
}

process DOWNLOAD_EGGNOG {
    tag 'eggnog'
    storeDir "${params.db_dir}/eggnog"
    output: path 'eggnog.db'
    script:
    """
    download_eggnog_data.py -y --data_dir . ${params.eggnog_download_args}
    """
    stub:
    "touch eggnog.db"
}
