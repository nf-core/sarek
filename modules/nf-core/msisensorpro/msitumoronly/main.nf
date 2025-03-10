process MSISENSORPRO_MSITUMORONLY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/msisensor-pro:1.2.0--hfc31af2_0' :
        'biocontainers/msisensor-pro:1.2.0--hfc31af2_0' }"

    input:
    tuple val(meta), path(tumor), path(tumor_index), path(intervals)
    path (fasta)
    path (msisensor_baseline)

    output:
    tuple val(meta), path("${prefix}")         , emit: output_report
    tuple val(meta), path("${prefix}_dis")     , emit: output_dis
    tuple val(meta), path("${prefix}_all")     , emit: output_all
    tuple val(meta), path("${prefix}_unstable"), emit: output_unstable
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def fasta = fasta ? "-g ${fasta}" : ""
    def intervals = intervals ? "-e ${intervals}" : ""
    """
    msisensor-pro \\
        pro \\
        -d ${msisensor_baseline} \\
        -t ${tumor} \\
        ${fasta} \\
        -o $prefix \\
        -b ${task.cpus} \\
        ${intervals} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msisensor-pro: \$(msisensor-pro 2>&1 | sed -nE 's/Version:\\sv([0-9]\\.[0-9])/\\1/ p')
    END_VERSIONS
    """
}
