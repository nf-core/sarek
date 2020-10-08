include { initOptions; saveFiles; getSoftwareName } from './../functions'

process GATK_GATHERBQSRREPORTS {
    label 'memory_singleCPU_2_task'
    label 'cpus_2'
    tag "${meta.id}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {
            if (it == "${meta.sample}.recal.table" && !params.skip_markduplicates) "preprocessing/${meta.sample}/markduplicates/${it}"
            else "preprocessing/${meta.sample}/mapped/${it}"
        }
 
    container "quay.io/biocontainers/gatk4-spark:4.1.8.1--0"

    conda (params.conda ? "bioconda::gatk4-spark=4.1.8.1" : null)

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