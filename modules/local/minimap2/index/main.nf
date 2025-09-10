process MINIMAP2_INDEX {
    tag "$fasta"
    label 'process_medium'

    conda "bioconda::minimap2=2.26"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.26--he4a0461_2' :
        'biocontainers/minimap2:2.26--he4a0461_2' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.mmi"), emit: index
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    minimap2 \\
        -x sr \\
        -t $task.cpus \\
        -d ${fasta}.mmi \\
        $args \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
