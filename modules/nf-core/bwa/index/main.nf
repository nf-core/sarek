process BWA_INDEX {
    tag "$fasta"
    label 'process_single'

    conda "bioconda::bwa=0.7.17"
    container { NfcoreTemplate.getContainer(
        docker: 'biocontainers/bwa:0.7.17--hed695b0_7',
        singularity: 'https://depot.galaxyproject.org/singularity/bwa:0.7.17--hed695b0_7',
        registry: 'quay.io',
        engine: workflow.containerEngine,
        use_full_uri: task.ext.container_full_uri ?: false,
    ) }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(bwa) , emit: index
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir bwa
    bwa \\
        index \\
        $args \\
        -p bwa/${fasta.baseName} \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir bwa

    touch bwa/genome.amb
    touch bwa/genome.ann
    touch bwa/genome.bwt
    touch bwa/genome.pac
    touch bwa/genome.sa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}
