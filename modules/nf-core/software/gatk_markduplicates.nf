process GATK_MARKDUPLICATES {
    label 'cpus_16'
    //tag "${meta.id}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {
            if (it == "${meta.sample}.bam.metrics") "Reports/${meta.sample}/MarkDuplicates/${it}"
            else "Preprocessing/${meta.sample}/DuplicatesMarked/${it}"
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

    if (params.no_gatk_spark)
    """
    gatk --java-options ${markdup_java_options} \
        MarkDuplicates \
        --MAX_RECORDS_IN_RAM 50000 \
        --INPUT ${meta.sample}.bam \
        --METRICS_FILE ${meta.sample}.bam.metrics \
        --TMP_DIR . \
        --ASSUME_SORT_ORDER coordinate \
        --CREATE_INDEX true \
        --OUTPUT ${meta.sample}.md.bam
    mv ${meta.sample}.md.bai ${meta.sample}.md.bam.bai
    """
    else
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
}