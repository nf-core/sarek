process GATK4SPARK_MARKDUPLICATES {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/49/498aea9c9bcaf736b9fb2a01366c1b7b38ccc0d38143178afc325d6a93241447/data'
        : 'community.wave.seqera.io/library/gatk4-spark:4.6.2.0--8b5cd67ee60a714e'}"

    input:
    tuple val(meta), path(bam)
    path fasta
    path fasta_fai
    path dict

    output:
    tuple val(meta), path("${prefix}"), emit: output
    tuple val(meta), path("${prefix}.bai"), emit: bam_index, optional: true
    tuple val(meta), path("*.metrics"), emit: metrics, optional: true
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.bam"
    def input_list = bam.collect { bam_ -> "--input ${bam_}" }.join(' ')

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK MarkDuplicatesSpark] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        MarkDuplicatesSpark \\
        ${input_list} \\
        --output ${prefix} \\
        --reference ${fasta} \\
        --spark-master local[${task.cpus}] \\
        --tmp-dir . \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.bam"
    """
    touch ${prefix}
    touch ${prefix}.bai
    touch ${prefix}.metrics
    """
}
