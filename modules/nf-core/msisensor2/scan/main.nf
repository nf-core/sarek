process MSISENSOR2_SCAN {
    tag "${fasta}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/msisensor2:0.1--hd03093a_0'
        : 'biocontainers/msisensor2:0.1--hd03093a_0'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.scan"), emit: scan
    path "versions.yml",             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def inputs = fasta.sort().collect { "-d ${it}" }.join(" ")
    """
    msisensor2 scan \\
        ${args} \\
        ${inputs} \\
        -o ${prefix}.scan

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msisensor2: \$(echo \$(msisensor2 2> >(grep Version) | sed 's/Version: v//g'))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.scan

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msisensor2: \$(echo \$(msisensor2 2> >(grep Version) | sed 's/Version: v//g'))
    END_VERSIONS
    """
}
