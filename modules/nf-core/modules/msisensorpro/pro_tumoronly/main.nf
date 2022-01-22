process MSISENSORPRO_MSI_SOMATIC {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda:: msisensor-pro=1.2.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/msisensor-pro:1.2.0--hfc31af2_0' :
        'quay.io/biocontainers/msisensor-pro:1.2.0--hfc31af2_0' }"

    input:
    tuple val(meta), path(input), path(index), val(metascan), path(homopolymers)

    output:
    tuple val(meta), path("${prefix}")         , emit: output
    tuple val(meta), path("${prefix}_dis")     , emit: output_dis
    tuple val(meta), path("${prefix}_germline"), emit: output_germline
    tuple val(meta), path("${prefix}_somatic") , emit: output_somatic
    path "versions.yml"                        , emit: versions

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    msisensor-pro \\
        msi \\
        -d $homopolymers \\
        -n $normal_bam \\
        -t $tumor_bam \\
        -o $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msisensor: \$(msisensor 2>&1 | sed -nE 's/Version:\\sv([0-9]\\.[0-9])/\\1/ p')
    END_VERSIONS
    """
}
