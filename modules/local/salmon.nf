//
// Salmon: decoy-aware selective-alignment quantification.
//
// Two index shapes:
//   - CDS-only: fast, small, always available. Default at `huge` genome size
//     (a decoy-aware index at 21.7 Gb is a comparably heavy build to the
//     HISAT2 index itself).
//   - decoy-aware: transcriptome + full genome as decoy, absorbs reads from
//     unannotated/repetitive regions that would otherwise be forced onto the
//     nearest transcript. Default everywhere the genome is small enough that
//     the extra build cost is cheap. See docs/decisions.md.
//
// Two transcript-input shapes:
//   - one haplotype's CDS: ordinary gene-level quantification.
//   - two haplotypes' CDS concatenated: Salmon's equivalence-class EM resolves
//     the ~99%-identical multi-mapping between alleles that HISAT2 + counting
//     cannot (see docs/haplotypes.md). IDs are checked disjoint before
//     concatenation -- the same "chr1a means two things" trap, applied to CDS.
//

process SALMON_CHECK_DISJOINT_IDS {
    label 'process_single'
    tag "${fastas.collect{it.name}.join(',')}"

    input:
    path fastas

    output:
    path 'combined.cds.fa', emit: fasta
    path 'versions.yml'   , emit: versions

    script:
    """
    check_disjoint_fasta_ids.py --out combined.cds.fa ${fastas}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    cat ${fastas} > combined.cds.fa
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}

process SALMON_INDEX {
    label 'process_high'
    tag "${transcript_fasta.name}"

    input:
    path transcript_fasta
    path genome_fasta     // NO_FILE if no decoy requested
    val  use_decoy

    output:
    path 'salmon_index', emit: index
    path 'versions.yml', emit: versions

    script:
    if( use_decoy && genome_fasta.name != 'NO_FILE' ) {
        """
        # Decoy-aware index: transcriptome sequences first, then the genome as
        # a decoy set, with a decoys.txt naming which sequence names are decoys
        # (the genome's own sequence names) -- the standard Salmon recipe.
        grep '^>' '${genome_fasta}' | sed 's/^>//' | cut -d' ' -f1 > decoys.txt
        cat '${transcript_fasta}' '${genome_fasta}' > gentrome.fa

        salmon index -t gentrome.fa -d decoys.txt -i salmon_index \\
            -k 31 -p ${task.cpus}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            salmon: \$(salmon --version | sed 's/salmon //')
        END_VERSIONS
        """
    } else {
        """
        salmon index -t '${transcript_fasta}' -i salmon_index -k 31 -p ${task.cpus}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            salmon: \$(salmon --version | sed 's/salmon //')
        END_VERSIONS
        """
    }

    stub:
    """
    mkdir -p salmon_index && touch salmon_index/info.json
    echo '"${task.process}": {salmon: stub}' > versions.yml
    """
}

process SALMON_QUANT {
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path("${meta.id}"), emit: quant_dir
    tuple val(meta), path("${meta.id}/quant.sf"), emit: quant_sf
    path 'versions.yml', emit: versions

    script:
    // -l A: Salmon's own automatic library-type detection. Independent of, and
    // reported alongside, this pipeline's junction-motif strandedness call
    // (see INFER_STRANDEDNESS) -- agreement between the two is a QC signal,
    // surfaced in MultiQC rather than silently trusted.
    def reads_arg = meta.single_end
        ? "-r ${reads[0]}"
        : "-1 ${reads[0..<(reads.size()/2 as int)].join(' ')} -2 ${reads[(reads.size()/2 as int)..<reads.size()].join(' ')}"
    """
    salmon quant -i '${index}' -l A ${reads_arg} \\
        -p ${task.cpus} --validateMappings --gcBias --seqBias \\
        -o ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version | sed 's/salmon //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}
    printf 'Name\\tLength\\tEffectiveLength\\tTPM\\tNumReads\\n' > ${meta.id}/quant.sf
    echo '"${task.process}": {salmon: stub}' > versions.yml
    """
}

process SALMON_TO_GENE_MATRIX {
    label 'process_low'
    tag 'salmon'

    input:
    val   quant_pairs   // list of "sample_id=path" strings
    path  quant_files    // the quant.sf files themselves, staged alongside
    path  gff3            // NO_FILE if no mapping available

    output:
    path 'counts_salmon.tsv', emit: counts
    path 'tpm_salmon.tsv'   , emit: tpm
    path 'versions.yml'     , emit: versions

    script:
    def gff_arg = gff3.name != 'NO_FILE' ? "--gff3 '${gff3}'" : ''
    """
    salmon_to_gene_matrix.py \\
        --quant ${quant_pairs.join(' ')} \\
        ${gff_arg} \\
        --out-counts counts_salmon.tsv --out-tpm tpm_salmon.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf 'gene_id\\n' > counts_salmon.tsv
    printf 'gene_id\\n' > tpm_salmon.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}

process SPLIT_HAPLOTYPE_COUNTS {
    label 'process_single'
    tag "${matrix.name}"

    input:
    path matrix
    path allele_pairs   // from haplotype_pairing's SYRI_TO_ALLELE_TABLE
    val  label           // 'counts' or 'tpm', for output naming only

    output:
    path "summed_${label}.tsv"          , emit: summed
    path "allele_specific_${label}.tsv" , emit: allele_specific
    path "unpaired_${label}.tsv"        , emit: unpaired
    path 'versions.yml'                 , emit: versions

    script:
    """
    split_haplotype_counts.py \\
        --matrix '${matrix}' --pairs '${allele_pairs}' \\
        --out-summed summed_${label}.tsv \\
        --out-allele-specific allele_specific_${label}.tsv \\
        --out-unpaired unpaired_${label}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch summed_${label}.tsv allele_specific_${label}.tsv unpaired_${label}.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
