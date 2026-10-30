process CNVKIT_CALL {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/cnvkit:0.9.12--pyhdfd78af_0'
        : 'quay.io/biocontainers/cnvkit:0.9.12--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(cns), path(vcf)

    output:
    tuple val(meta), path("*.cns"), emit: cns
    tuple val("${task.process}"), val('cnvkit'), eval('cnvkit.py version | sed -e "s/cnvkit v//g"'), emit: versions_cnvkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcf_cmd = vcf ? "-v ${vcf}" : ""
    """
    cnvkit.py call \\
        ${cns} \\
        ${vcf_cmd} \\
        ${args} \\
        -o ${prefix}.cns
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cns
    """
}
