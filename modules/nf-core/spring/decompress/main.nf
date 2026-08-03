process SPRING_DECOMPRESS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f6/f67f27c8cb2d1a149564f1a10f5f2b7a6acfa87ef3d3d27d2d8752dbe95e6acf/data' :
        'community.wave.seqera.io/library/spring:1.1.1--911a17b4ccfb85ee' }"

    input:
    tuple val(meta), path(spring)
    val(write_one_fastq_gz)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    tuple val("${task.process}"), val('spring'), val('1.1.1'), topic: versions, emit: versions_spring
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = write_one_fastq_gz ? "-o ${prefix}.fastq.gz" : "-o ${prefix}_R1.fastq.gz ${prefix}_R2.fastq.gz"
    """
    spring \\
        -d \\
        -g \\
        -t ${task.cpus} \\
        $args \\
        -i ${spring} \\
        ${output}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = write_one_fastq_gz ? "echo '' | gzip > ${prefix}.fastq.gz" : "echo '' | gzip > ${prefix}_R1.fastq.gz; echo '' | gzip > ${prefix}_R2.fastq.gz"
    """
    ${output}
    """
}
