process TABIX_TABIX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/htslib:1.19.1--h81da01d_1' :
        'biocontainers/htslib:1.19.1--h81da01d_1' }"

    input:
    tuple val(meta), path(tab)

    output:
    tuple val(meta), path("*.tbi"), optional:true, emit: tbi
    tuple val(meta), path("*.csi"), optional:true, emit: csi
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    tabix $args $tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${tab}.tbi
    touch ${tab}.csi
    cat <<-END_VERSIONS > versions.yml

    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
