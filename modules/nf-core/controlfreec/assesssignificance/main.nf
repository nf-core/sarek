process CONTROLFREEC_ASSESSSIGNIFICANCE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/control-freec:11.6--h1b792b2_1' :
        'biocontainers/control-freec:11.6--h1b792b2_1' }"

    input:
    tuple val(meta), path(cnvs), path(ratio)

    output:
    tuple val(meta), path("*.p.value.txt"), emit: p_value_txt
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '11.6b' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    cat \$(which assess_significance.R) | R --slave --args ${cnvs} ${ratio}

    mv *.p.value.txt ${prefix}.p.value.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        controlfreec: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '11.6b' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.p.value.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        controlfreec: $VERSION
    END_VERSIONS
    """
}
