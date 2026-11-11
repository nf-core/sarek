process GAWK {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a1/a125c778baf3865331101a104b60d249ee15fe1dca13bdafd888926cc5490a34/data' :
        'community.wave.seqera.io/library/gawk:5.3.1--e09efb5dfc4b8156' }"

    input:
    tuple val(meta), path(input, arity: '0..*')
    path(program_file)
    val(disable_redirect_output)

    output:
    tuple val(meta), path("*.${suffix}"), emit: output
    tuple val("${task.process}"), val('gawk'), eval("awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//'"), topic: versions, emit: versions_gawk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args  ?: '' // args is used for the main arguments of the tool
    def args2 = task.ext.args2 ?: '' // args2 is used to specify a program when no program file has been given
    prefix    = task.ext.prefix ?: "${meta.id}"
    suffix    = task.ext.suffix ?: "${input.collect{ file -> file.getExtension()}.get(0)}" // use the first extension of the input files

    program    = program_file ? "-f ${program_file}" : "${args2}"
    lst_gz     = input.findResults{ file -> file.getExtension().endsWith("gz") ? file.toString() : null }
    unzip      = lst_gz ? "gunzip -q -f ${lst_gz.join(" ")}" : ""
    input_cmd  = input.collect { file -> file.toString() - ~/\.gz$/ }.join(" ")
    output_cmd = suffix.endsWith("gz") ? "| gzip > ${prefix}.${suffix}" : "> ${prefix}.${suffix}"
    output     = disable_redirect_output ? "" : output_cmd
    cleanup    = lst_gz ? "rm ${lst_gz.collect{ file -> file - ~/\.gz$/ }.join(" ")}" : ""

    input.collect{ file ->
        assert file.name != "${prefix}.${suffix}" : "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }

    """
    ${unzip}

    awk \\
        ${args} \\
        ${program} \\
        ${input_cmd} \\
        ${output}

    ${cleanup}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = task.ext.suffix ?: "${input.collect{ file -> file.getExtension()}.get(0)}"
    def create_cmd = suffix.endsWith("gz") ? "echo '' | gzip >" : "touch"

    """
    ${create_cmd} ${prefix}.${suffix}
    """
}
