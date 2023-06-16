process DRAGMAP_HASHTABLE {
    tag "$fasta"
    label 'process_high'

    conda "bioconda::dragmap=1.2.1"
    container { NfcoreTemplate.getContainer(
        docker: 'biocontainers/dragmap:1.2.1--h72d16da_1',
        singularity: 'https://depot.galaxyproject.org/singularity/dragmap:1.2.1--h72d16da_1',
        registry: 'quay.io',
        use_full_uri: task.ext.container_full_uri ?: false
    )}

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("dragmap")    , emit: hashmap
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir dragmap
    dragen-os \\
        --build-hash-table true \\
        --ht-reference $fasta \\
        --output-directory dragmap \\
        $args \\
        --ht-num-threads $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragmap: \$(echo \$(dragen-os --version 2>&1))
    END_VERSIONS
    """
}
