include { initOptions; saveFiles; getSoftwareName } from './../functions'

params.options = [:]
def options    = initOptions(params.options)

process GATK_MARKDUPLICATES {
    label 'cpus_16'

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
        tuple val(meta), path("${meta.sample}.bam"), path("${meta.sample}.bam.bai")

    output:
        tuple val(meta), path("${meta.sample}.md.bam"), path("${meta.sample}.md.bam.bai"), emit: bam
        val meta,                                                                          emit: tsv
        path "${meta.sample}.bam.metrics", optional : true,                                emit: report
          
    script:
    markdup_java_options = task.memory.toGiga() > 8 ? params.markdup_java_options : "\"-Xms" +  (task.memory.toGiga() / 2).trunc() + "g -Xmx" + (task.memory.toGiga() - 1) + "g\""
    metrics = 'markduplicates' in params.skip_qc ? '' : "-M ${meta.sample}.bam.metrics"

    if (params.use_gatk_spark)
    """
    gatk --java-options ${markdup_java_options} \
        MarkDuplicatesSpark \
        -I ${meta.sample}.bam \
        -O ${meta.sample}.md.bam \
        ${metrics} \
        --tmp-dir . \
        --create-output-bam-index true \
        --spark-master local[${task.cpus}]
    """
    else
    """
    gatk --java-options ${markdup_java_options} \
        MarkDuplicates \
        --INPUT ${meta.sample}.bam \
        --METRICS_FILE ${meta.sample}.bam.metrics \
        --TMP_DIR . \
        --ASSUME_SORT_ORDER coordinate \
        --CREATE_INDEX true \
        --OUTPUT ${meta.sample}.md.bam
    mv ${meta.sample}.md.bai ${meta.sample}.md.bam.bai
    """
}