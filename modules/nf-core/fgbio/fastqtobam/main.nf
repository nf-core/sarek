process FGBIO_FASTQTOBAM {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    // --version argument gives the wrong version
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:2.1.0--hdfd78af_0' :
        'biocontainers/fgbio:2.1.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam") , emit: bam , optional: true
    tuple val(meta), path("*.cram"), emit: cram, optional: true
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "bam"
    def sample_name = args.contains("--sample") ? "" : "--sample ${prefix}"
    def library_name = args.contains("--library") ? "" : "--library ${prefix}"

    def VERSION = '2.1.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """

    fgbio \\
        --tmp-dir=. \\
        FastqToBam \\
        ${args} \\
        --input ${reads} \\
        --output ${prefix}.${suffix} \\
        ${sample_name} \\
        ${library_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "bam"

    def VERSION = '2.1.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: $VERSION
    END_VERSIONS
    """
}
