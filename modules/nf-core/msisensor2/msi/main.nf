process MSISENSOR2_MSI {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/msisensor2:0.1--hd03093a_0'
        : 'quay.io/biocontainers/msisensor2:0.1--hd03093a_0'}"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bam_index)
    tuple val(meta2), path(models)

    output:
    tuple val(meta), path("${prefix}"), emit: msi
    tuple val(meta), path("${prefix}_dis"), emit: distribution
    tuple val(meta), path("${prefix}_somatic"), emit: somatic
    tuple val("${task.process}"), val('msisensor2'), eval("msisensor2 2> >(grep Version) | sed 's/Version: v//g'"), topic: versions, emit: versions_msisensor2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    msisensor2 msi \\
        -b ${task.cpus} \\
        ${args} \\
        -M ${models} \\
        -t ${tumor_bam} \\
        -o ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}
    touch ${prefix}_dis
    touch ${prefix}_somatic
    """
}
