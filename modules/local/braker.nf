//
// BRAKER3: RNA-seq guided gene prediction
// Input: genome FASTA + RNA-seq BAM (from HISAT2)
// Output: GFF3 gene models (braker.gff3), protein FASTA (braker.aa), coding sequences (braker.cds)
//
process BRAKER3 {
    label 'process_high'
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bam_index)
    path genome_fasta
    path genome_fai
    val strategy

    output:
    tuple val(meta), path('braker.gff3'), path('braker.aa'), path('braker.cds'), emit: gff
    path 'versions.yml', emit: versions

    script:
    def max_intron = params.braker_max_intron ?: 100000
    def min_intron = params.braker_min_intron ?: 20
    def threads = task.cpus
    def memory = "${params.braker_memory ?: '16.GB'}"

    // Prepare proteins if provided
    def proteins = params.braker_proteins ?: ''
    def proteins_arg = proteins ? "--proteins=${proteins}" : ''

    // Softmasking already done in GENOME_PREP, so we use --softmasking
    def softmasking = params.braker_softmasking ?: true
    def softmask_arg = softmasking ? '--softmasking' : ''

    // AUGUSTUS species hint
    def augustus_species = params.braker_augustus_species ?: ''
    def augustus_arg = augustus_species ? "--species=${augustus_species}" : ''

    // UTR training
    def utr = params.braker_utr ?: false
    def utr_arg = utr ? '--UTR=on' : ''

    // AUGUSTUS config path
    def augustus_config_path = params.augustus_config_path ?: ''
    def augustus_config_arg = augustus_config_path ? "--AUGUSTUS_CONFIG_PATH=${augustus_config_path}" : ''

    // BAM input handling - BRAKER3 expects sorted BAM
    def bam_input = bam

    """
    set -euo pipefail

    # Setup working directory
    mkdir -p braker_out
    cd braker_out

    # Link input files
    ln -sf ${genome_fasta} genome.fa
    ln -sf ${genome_fai} genome.fa.fai
    ln -sf ${bam_input} rnaseq.bam
    ln -sf ${bam_index} rnaseq.bam.bai

    # Run BRAKER3
    braker.pl \\
        --genome=genome.fa \\
        --bam=rnaseq.bam \\
        --threads=${threads} \\
        --maxIntronLen=${max_intron} \\
        --minIntronLen=${min_intron} \\
        ${softmask_arg} \\
        ${augustus_arg} \\
        ${utr_arg} \\
        ${proteins_arg} \\
        --workingdir=. \\
        --cores=${threads} \\
        --AUGUSTUS_CONFIG_PATH=/opt/conda/envs/braker/config/ \\
        2>&1 | tee braker.log

    # Copy outputs
    cp braker/braker.gtf ../braker.gtf
    cp braker/braker.aa ../braker.aa
    cp braker/braker.cds ../braker.cds

    # Convert GTF to GFF3
    gffread braker/braker.gtf -T -o braker.gff3
    cp braker.gff3 ../braker.gff3

    # Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        braker: \$(braker.pl --version 2>&1 | head -1 | sed 's/.* //')
        augustus: \$(augustus --version 2>&1 | head -1 | sed 's/.* //')
        genemark: \$(gmhmme3 -v 2>&1 | head -1 | sed 's/.* //')
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS

    """

    stub:
    """
    mkdir -p braker_out
    echo '##gff-version 3' > braker.gff3
    echo 'chr1	BRAKER3	gene	100	500	.	+	.	ID=gene1' >> braker.gff3
    echo '>protein1' > braker.aa
    echo 'MKT...' >> braker.aa
    echo 'ATG...' > braker.cds
    echo '"${task.process}": {braker: stub, augustus: stub, genemark: stub}' > versions.yml
    """
}

process BRAKER3_MERGE {
    label 'process_low'
    tag "merge"

    input:
    path gff_files

    output:
    path 'braker_merged.gff3', emit: gff
    path 'versions.yml', emit: versions

    script:
    """
    cat ${gff_files} | grep -v '^##gff-version' | sort -k1,1 -k4,4n > braker_merged.gff3
    echo '##gff-version 3' | cat - braker_merged.gff3 > braker_merged.gff3.tmp && mv braker_merged.gff3.tmp braker_merged.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    echo '##gff-version 3' > braker_merged.gff3
    echo '"${task.process}": {gffread: stub}' > versions.yml
    """
}