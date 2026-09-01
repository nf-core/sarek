process CNVKIT_EXPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/cnvkit:0.9.12--pyhdfd78af_0'
        : 'quay.io/biocontainers/cnvkit:0.9.12--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(cns)

    output:
    tuple val(meta), path("${prefix}.${suffix}"), emit: output
    tuple val("${task.process}"), val('cnvkit'), eval('cnvkit.py version | sed -e "s/cnvkit v//g"'), emit: versions_cnvkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = task.ext.args.tokenize(" ")[0]
    """
    cnvkit.py export \\
        ${args} \\
        ${cns} \\
        -o ${prefix}.${suffix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = task.ext.args.tokenize(" ")[0]
    """
    touch ${prefix}.${suffix}
    """
}
