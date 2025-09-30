process YTE {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f9/f9c1bae08d5c03036185fde01db24f11173ab0d746887a3313ec00851b2fb5ae/data':
        'community.wave.seqera.io/library/pip_yte:93491093a59d72ba' }"

    input:
    tuple val(meta), path(template)
    path(map_file)
    val(map)

    output:
    tuple val(meta), path("*.yaml"), emit: rendered
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Use map_file if provided, otherwise use map to create key=value pairs for mapping command
    def mapping_cmd = map_file ? "--variable-file ${map_file}" : "--variables " + map.collect { k, v -> "${k}=${v}" }.join(' ')
    VERSION = "1.9.0" // WARN: Version information not provided by tool on CLI. Please update this string when bumping
    """
    yte ${mapping_cmd} < ${template} > ${prefix}.yaml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yte: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    VERSION = "1.9.0" // WARN: Version information not provided by tool on CLI. Please update this string when bumping
    """
    touch ${prefix}.yaml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yte: $VERSION
    END_VERSIONS
    """
}
