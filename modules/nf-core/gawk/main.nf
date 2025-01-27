process GAWK {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(input, arity: '0..*')
    path(program_file)

    output:
    tuple val(meta), path("${prefix}.${suffix}"), emit: output
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args  ?: '' // args is used for the main arguments of the tool
    def args2 = task.ext.args2 ?: '' // args2 is used to specify a program when no program file has been given
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = task.ext.suffix ?: "${input.collect{ it.getExtension()}.get(0)}" // use the first extension of the input files

    program   = program_file ? "-f ${program_file}" : "${args2}"
    lst_gz    = input.collect{ it.getExtension().endsWith("gz") }
    unzip     = lst_gz.contains(false) ? "" : "find ${input} -exec zcat {} \\; | \\"
    input_cmd = unzip ? "" : "${input}"

    """
    ${unzip}
    awk \\
        ${args} \\
        ${program} \\
        ${input_cmd} \\
        > ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = task.ext.suffix ?: "${input.getExtension()}"
    def create_cmd = suffix.endsWith("gz") ? "echo '' | gzip >" : "touch"

    """
    ${create_cmd} ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
