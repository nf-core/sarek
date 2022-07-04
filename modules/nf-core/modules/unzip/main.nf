process UNZIP {
    tag "$archive"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::p7zip=15.09" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/p7zip:15.09--h2d50403_4' :
        'quay.io/biocontainers/p7zip:15.09--h2d50403_4' }"

    input:
    tuple val(meta), path(archive)

    output:
    tuple val(meta), path("${archive.baseName}/"), emit: unzipped_archive
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    if ( archive instanceof List && archive.name.size > 1 ) { exit 1, "[UNZIP] error: 7za only accepts a single archive as input. Please check module input." }
    """
    7za \\
        e \\
        -o"${archive.baseName}"/ \\
        $args \\
        $archive

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        7za: \$(echo \$(7za --help) | sed 's/.*p7zip Version //; s/(.*//')
    END_VERSIONS
    """
}
