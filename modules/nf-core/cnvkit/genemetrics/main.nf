process CNVKIT_GENEMETRICS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/cnvkit:0.9.12--pyhdfd78af_0'
        : 'quay.io/biocontainers/cnvkit:0.9.12--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(cnr), path(cns)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('cnvkit'), eval('cnvkit.py version | sed -e "s/cnvkit v//g"'), emit: versions_cnvkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def segments = cns ? "--segment ${cns}" : ""

    """
    cnvkit.py \\
        genemetrics \\
        ${cnr} \\
        ${segments} \\
        --output ${prefix}.tsv \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.tsv
    """
}
