process CNVKIT_ANTITARGET {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::cnvkit=0.9.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvkit:0.9.9--pyhdfd78af_0':
        'quay.io/biocontainers/cnvkit:0.9.9--pyhdfd78af_0' }"

    input:
    path  targets

    output:
    path("*.bed")                 , emit: BED
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "antitarget"

    """
    cnvkit.py \\
        antitarget \\
        $targets \\
        --output ${prefix}.bed \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$(cnvkit.py version | sed -e "s/cnvkit v//g")
    END_VERSIONS
    """
}
