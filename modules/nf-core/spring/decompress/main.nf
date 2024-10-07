process SPRING_DECOMPRESS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/spring:1.1.1--h4ac6f70_2' :
        'biocontainers/spring:1.1.1--h4ac6f70_2' }"

    input:
    tuple val(meta), path(spring)
    val(write_one_fastq_gz)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def output = write_one_fastq_gz ? "-o ${prefix}.fastq.gz" : "-o ${prefix}_R1.fastq.gz ${prefix}_R2.fastq.gz"

    """
    spring \\
        -d \\
        -g \\
        -t ${task.cpus} \\
        $args \\
        -i ${spring} \\
        ${output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spring: ${VERSION}
    END_VERSIONS
    """
}
