// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BEDTOOLS_SORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52d3fbb86c48b8c39a989733ca97bfd84cf665a9ab279be7cfbbb36f5a53c46a/data'
        : 'community.wave.seqera.io/library/bedtools:2.31.1--897920cd303d6413' }"

    input:
    tuple val(meta), path(intervals)
    path genome_file
    val suffix

    output:
    tuple val(meta), path("*.${extension}"), emit: sorted
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed -e 's/bedtools v//g'"), topic: versions, emit: versions_bedtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genome_cmd = genome_file ? "-g ${genome_file}" : ""
    extension = suffix ?: intervals.extension
    if ("${intervals}" == "${prefix}.${extension}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    bedtools \\
        sort \\
        -i ${intervals} \\
        ${genome_cmd} \\
        ${args} \\
        > ${prefix}.${extension}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    extension = suffix ?: intervals.extension
    """
    touch ${prefix}.${extension}
    """
}
