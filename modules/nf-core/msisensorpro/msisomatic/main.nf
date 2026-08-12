process MSISENSORPRO_MSISOMATIC {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/msisensor-pro%3A1.3.0--hfef96ef_0':
            'biocontainers/msisensor-pro:1.3.0--hfef96ef_0' }"

    input:
    tuple val(meta), path(normal), path(normal_index), path(tumor), path(tumor_index), path(intervals)
    tuple val(meta2), path(fasta)
    path(msisensor_scan)

    output:
    tuple val(meta), path("${prefix}")         , emit: output_report
    tuple val(meta), path("${prefix}_dis")     , emit: output_dis
    tuple val(meta), path("${prefix}_germline"), emit: output_germline, optional: true
    tuple val(meta), path("${prefix}_somatic") , emit: output_somatic,  optional: true
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def fasta_cmd     = fasta     ? "-g ${fasta}"       : ""
    def intervals_cmd = intervals ? " -e ${intervals} " : ""

    """
    msisensor-pro \\
        msi \\
        -d ${msisensor_scan} \\
        -n ${normal} \\
        -t ${tumor} \\
        ${fasta_cmd} \\
        -o ${prefix} \\
        -b ${task.cpus} \\
        ${intervals_cmd} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msisensor-pro: \$(msisensor-pro 2>&1 | sed -nE 's/Version:\\sv([0-9]\\.[0-9])/\\1/ p')
    END_VERSIONS
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}
    touch ${prefix}_dis
    touch ${prefix}_germline
    touch ${prefix}_somatic

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msisensor-pro: \$(msisensor-pro 2>&1 | sed -nE 's/Version:\\sv([0-9]\\.[0-9])/\\1/ p')
    END_VERSIONS
    """
}
