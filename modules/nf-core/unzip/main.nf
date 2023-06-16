process UNZIP {
    tag "$archive"
    label 'process_single'

    conda "conda-forge::p7zip=16.02"
    container { NfcoreTemplate.getContainer(
        docker: 'biocontainers/p7zip:16.02',
        singularity: 'https://depot.galaxyproject.org/singularity/p7zip:16.02',
        registry: 'quay.io',
        use_full_uri: task.ext.container_full_uri ?: false
    )}

    input:
    tuple val(meta), path(archive)

    output:
    tuple val(meta), path("${prefix}/"), emit: unzipped_archive
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    if ( archive instanceof List && archive.name.size > 1 ) { exit 1, "[UNZIP] error: 7za only accepts a single archive as input. Please check module input." }

    prefix = task.ext.prefix ?: ( meta.id ? "${meta.id}" : archive.baseName)
    """
    7za \\
        x \\
        -o"${prefix}"/ \\
        $args \\
        $archive

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        7za: \$(echo \$(7za --help) | sed 's/.*p7zip Version //; s/(.*//')
    END_VERSIONS
    """
}
