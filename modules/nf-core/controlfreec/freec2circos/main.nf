process CONTROLFREEC_FREEC2CIRCOS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/control-freec:11.6b--hde5307d_3'
        : 'quay.io/biocontainers/control-freec:11.6b--hde5307d_3'}"

    input:
    tuple val(meta), path(ratio)

    output:
    tuple val(meta), path("*.circos.txt"), emit: circos
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('controlfreec'), val("11.6b"), emit: versions_controlfreec, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    freec2circos.pl -f ${ratio} ${args} > ${prefix}.circos.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.circos.txt
    """
}
