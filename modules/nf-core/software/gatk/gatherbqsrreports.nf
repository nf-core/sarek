include { initOptions; saveFiles; getSoftwareName } from './../functions'

params.options = [:]
def options    = initOptions(params.options)

process GATK_GATHERBQSRREPORTS {
    label 'memory_singleCPU_2_task'
    label 'cpus_2'

    tag "${meta.id}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }
 
    conda (params.enable_conda ? "bioconda::gatk4-spark=4.1.8.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4-spark:4.1.8.1--0"
    } else {
        container "quay.io/biocontainers/gatk4-spark:4.1.8.1--0"
    }

    input:
        tuple val(meta), path(recal)

    output:
        tuple val(meta), path("${meta.sample}.recal.table"), emit: table
        path "${meta.sample}.recal.table",                   emit: report
        val meta,                                            emit: tsv
        
    script:
    input = recal.collect{"-I ${it}"}.join(' ')
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        GatherBQSRReports \
        ${input} \
        -O ${meta.sample}.recal.table \
    """
}