process MSISENSOR_MSI {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::msisensor=0.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/msisensor:0.5--hb3646a4_2' :
        'quay.io/biocontainers/msisensor:0.5--hb3646a4_2' }"

    input:
    tuple val(meta), path(normal_bam), path(normal_bai), path(tumor_bam), path(tumor_bai), val(metascan), path(homopolymers)

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
    msisensor \\
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
