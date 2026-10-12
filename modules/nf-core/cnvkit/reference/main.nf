process CNVKIT_REFERENCE {
    tag "${fasta}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/cnvkit:0.9.12--pyhdfd78af_0'
        : 'quay.io/biocontainers/cnvkit:0.9.12--pyhdfd78af_0'}"

    input:
    path fasta
    path targets
    path antitargets

    output:
    path "*.cnn",        emit: cnn
    tuple val("${task.process}"), val('cnvkit'), eval('cnvkit.py version | sed -e "s/cnvkit v//g"'), emit: versions_cnvkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: targets.BaseName

    """
    cnvkit.py \\
        reference \\
        --fasta ${fasta} \\
        --targets ${targets} \\
        --antitargets ${antitargets} \\
        --output ${prefix}.reference.cnn \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: targets.BaseName

    """
    touch ${prefix}.reference.cnn
    """
}
