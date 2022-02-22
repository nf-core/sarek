process MSISENSORPRO_SCAN {
    tag "$fasta"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::msisensor-pro=1.1.a" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/msisensor-pro:1.1.a--hb3646a4_0' :
        'quay.io/biocontainers/msisensor-pro:1.1.a--hb3646a4_0' }"

    input:
    path  fasta

    output:
    tuple val(meta), path("*.list"), emit: list
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args  ?: ''
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    """
    msisensor-pro scan \\
        -d $fasta \\
        -o ${fasta.baseName}.list \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msisensor-pro: \$(msisensor-pro 2>&1 | sed -nE 's/Version:\\sv([0-9]\\.[0-9])/\\1/ p')
    END_VERSIONS
    """
}
