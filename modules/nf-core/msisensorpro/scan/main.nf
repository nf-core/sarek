process MSISENSORPRO_SCAN {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::msisensor-pro=1.2.0"
    container { NfcoreTemplate.getContainer(
        docker: 'biocontainers/msisensor-pro:1.2.0--hfc31af2_0',
        singularity: 'https://depot.galaxyproject.org/singularity/msisensor-pro:1.2.0--hfc31af2_0',
        registry: 'quay.io',
        use_full_uri: task.ext.container_full_uri ?: false
    )}

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.list"), emit: list
    path "versions.yml"            , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msisensor-pro: \$(msisensor-pro 2>&1 | sed -nE 's/Version:\\sv([0-9]\\.[0-9])/\\1/ p')
    END_VERSIONS
    """
}
