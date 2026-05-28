process CONTROLFREEC_ASSESSSIGNIFICANCE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/control-freec:11.6b--hde5307d_3'
        : 'quay.io/biocontainers/control-freec:11.6b--hde5307d_3'}"

    input:
    tuple val(meta), path(cnvs), path(ratio)

    output:
    tuple val(meta), path("*.p.value.txt"), emit: p_value_txt
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('controlfreec'), val("11.6b"), emit: versions_controlfreec, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat \$(which assess_significance.R) | R --slave --args ${cnvs} ${ratio}

    mv *.p.value.txt ${prefix}.p.value.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.p.value.txt
    """
}
