process DRAGMAP_HASHTABLE {
    tag "$fasta"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::dragmap=1.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dragmap:1.3.0--h72d16da_1':
        'quay.io/biocontainers/dragmap:1.3.0--h72d16da_1' }"

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
