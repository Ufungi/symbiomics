process INSPECT_VARIATION_FILE {
    label 'process_single'
    tag "${variation_file.name}"

    input:
    path variation_file
    val  accession

    output:
    path 'resequencing.vcf', emit: vcf, optional: true
    path 'versions.yml'    , emit: versions

    script:
    def acc = accession ? "--accession '${accession}'" : ''
    """
    inspect_variation_file.py --variation-file '${variation_file}' ${acc} --out resequencing.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch resequencing.vcf
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
