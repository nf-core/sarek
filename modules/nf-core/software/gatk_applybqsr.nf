process GATK_APPLYBQSR {
    label 'memory_singleCPU_2_task'
    label 'cpus_2'

    tag "${meta.id}-${interval.baseName}"

    input:
        tuple val(meta), path(bam), path(bai), path(recalibrationReport), file(interval)
        path dict
        path fasta
        path fai

    output:
        tuple val(meta), path("${prefix}${meta.sample}.recal.bam") 

    script:
    prefix = params.no_intervals ? "" : "${interval.baseName}_"
    options_intervals = params.no_intervals ? "" : "-L ${interval}"
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        ApplyBQSR \
        -R ${fasta} \
        --input ${bam} \
        --output ${prefix}${meta.sample}.recal.bam \
        ${options_intervals} \
        --bqsr-recal-file ${recalibrationReport}
    """
}