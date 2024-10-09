process CNVKIT_REFERENCE {
    tag "$fasta"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvkit:0.9.11--pyhdfd78af_0':
        'biocontainers/cnvkit:0.9.11--pyhdfd78af_0' }"

    input:
    path fasta
    path targets
    path antitargets

    output:
    path "*.cnn"       , emit: cnn
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: targets.BaseName

    """
    cnvkit.py \\
        reference \\
        --fasta $fasta \\
        --targets $targets \\
        --antitargets $antitargets \\
        --output ${prefix}.reference.cnn \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$(cnvkit.py version | sed -e "s/cnvkit v//g")
    END_VERSIONS
    """
}
