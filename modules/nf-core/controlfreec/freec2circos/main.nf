process CONTROLFREEC_FREEC2CIRCOS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/control-freec:11.6b--hdbdd923_0' :
        'biocontainers/control-freec:11.6b--hdbdd923_0' }"

    input:
    tuple val(meta), path(ratio)

    output:
    tuple val(meta), path("*.circos.txt"), emit: circos
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '11.6b' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    freec2circos.pl -f ${ratio} ${args} > ${prefix}.circos.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        controlfreec: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '11.6b' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.circos.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        controlfreec: $VERSION
    END_VERSIONS
    """
}
