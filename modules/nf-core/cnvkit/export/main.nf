process CNVKIT_EXPORT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvkit:0.9.10--pyhdfd78af_0':
        'biocontainers/cnvkit:0.9.10--pyhdfd78af_0' }"

    input:
    tuple val(meta) , path(cns)

    output:
    tuple val(meta), path("${prefix}.${suffix}"), emit: output
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = task.ext.args.tokenize(" ")[0]
    """
    cnvkit.py export \\
        $args \\
        $cns \\
        -o ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
    END_VERSIONS
    """
}
