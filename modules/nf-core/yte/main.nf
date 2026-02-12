process YTE {
    tag "${meta.id}"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a4/a427620cd214cbc7f48c8d43fa818be59c6c4ea4a622331d96a8b063a05335b0/data'
        : 'community.wave.seqera.io/library/yte:1.9.4--2a362f82cd32b54a'}"

    input:
    tuple val(meta), path(template), path(map_file), val(map)

    output:
    tuple val(meta), path("*.yaml"), emit: rendered
    tuple val("${task.process}"), val('yte'), eval("echo $VERSION"), topic: versions, emit: versions_yte

    when:
    task.ext.when == null || task.ext.when

    script:
    // No args because tool does not accept args, only stdin/stdout
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Use map_file if provided, otherwise use map to create key=value pairs for mapping command
    def mapping_cmd = map_file ? "--variable-file ${map_file}" : "--variables " + map.collect { k, v -> "${k}=${v}" }.join(' ')
    VERSION = "1.9.4"
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping
    """
    yte ${mapping_cmd} < ${template} > ${prefix}.yaml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    VERSION = "1.9.4"
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping
    """
    touch ${prefix}.yaml
    """
}
