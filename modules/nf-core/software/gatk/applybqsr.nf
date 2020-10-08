include { initOptions; saveFiles; getSoftwareName } from './../functions'

process GATK_APPLYBQSR {
    label 'memory_singleCPU_2_task'
    label 'cpus_2'

    tag "${meta.id}-${interval.baseName}"

    container "quay.io/biocontainers/gatk4-spark:4.1.8.1--0"

    conda (params.conda ? "bioconda::gatk4-spark=4.1.8.1" : null)

    input:
        tuple val(meta), path(bam), path(bai), path(recalibrationReport), path(interval)
        path dict
        path fasta
        path fai

    output:
        tuple val(meta), path("${prefix}${meta.sample}.recal.bam") , emit: bam
        val meta,                                                    emit: tsv


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