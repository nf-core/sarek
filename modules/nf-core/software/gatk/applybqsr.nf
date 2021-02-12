include { initOptions; saveFiles; getSoftwareName } from './../functions'

params.options = [:]
def options    = initOptions(params.options)

process GATK_APPLYBQSR {
    label 'memory_singleCPU_2_task'
    label 'cpus_2'

    tag "${meta.id}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::gatk4=4.1.9.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.1.9.0--py39_0"
    } else {
        container "quay.io/biocontainers/gatk4:4.1.9.0--py39_0"
    }

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