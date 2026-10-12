process MSISENSORPRO_SCAN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/msisensor-pro%3A1.3.0--hfef96ef_0':
            'quay.io/biocontainers/msisensor-pro:1.3.0--hfef96ef_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.list"), emit: list
    tuple val("${task.process}"), val('msisensor-pro'), eval("msisensor-pro --version 2>&1 | sed -nE 's/Version:\\s*v//p'") , emit: versions_msisensorpro, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    msisensor-pro \\
        scan \\
        -d $fasta \\
        -o ${prefix}.msisensor_scan.list \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.msisensor_scan.list
    """
}
