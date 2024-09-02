/**
 * The aim of this process is to re-index the bam file without the duplicate, supplementary, unmapped etc, for goleft/indexcov
 * It creates a BAM containing only a header (so indexcov can get the sample name) and a BAM index were low quality reads, supplementary etc, have been removed
 */
process SAMTOOLS_REINDEX_BAM {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0' :
        'biocontainers/samtools:1.20--h50ea8bc_0' }"

    input:
    tuple val(meta), path(input), path(input_index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("${meta.id}.reindex.bam"), path("${meta.id}.reindex.bam.bai"),emit: output
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def reference = fasta ? "--reference ${fasta}" : ""
    """
    # write header only
    samtools \\
        view \\
        --header-only \\
        --threads ${task.cpus} \\
        -O BAM \\
        -o "${meta.id}.reindex.bam" \\
        ${reference} \\
        ${input}

    # write BAM index only, remove unmapped, supplementary, etc...
    samtools \\
        view \\
        --uncompressed \\
        --write-index \\
        --threads ${task.cpus} \\
        -O BAM \\
        -o "/dev/null##idx##${meta.id}.reindex.bam.bai" \\
        ${reference} \\
        ${args} \\
        ${input}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
