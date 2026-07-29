process GATK4SPARK_APPLYBQSR {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/49/498aea9c9bcaf736b9fb2a01366c1b7b38ccc0d38143178afc325d6a93241447/data'
        : 'community.wave.seqera.io/library/gatk4-spark:4.6.2.0--8b5cd67ee60a714e'}"

    input:
    tuple val(meta), path(input), path(input_index), path(bqsr_table), path(intervals)
    path fasta
    path fai
    path dict

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam, optional: true
    tuple val(meta), path("${prefix}*bai"), emit: bai, optional: true
    tuple val(meta), path("${prefix}.cram"), emit: cram, optional: true
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // suffix can only be bam or cram, cram being the sensible default
    def suffix = task.ext.suffix && task.ext.suffix == "bam" ? "bam" : "cram"
    def interval_command = intervals ? "--intervals ${intervals}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK ApplyBQSRSpark] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        ApplyBQSRSpark \\
        --input ${input} \\
        --output ${prefix}.${suffix} \\
        --reference ${fasta} \\
        --bqsr-recal-file ${bqsr_table} \\
        ${interval_command} \\
        --spark-master local[${task.cpus}] \\
        --tmp-dir . \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "cram"
    """
    touch ${prefix}.${suffix}
    if [[ ${suffix} == bam ]]; then
        touch ${prefix}.${suffix}.bai
    fi
    """
}
