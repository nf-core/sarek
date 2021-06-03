/*
========================================================================================
    MARKDUPLICATES
========================================================================================
*/

params.markduplicates_options = [:]

include { GATK4_MARKDUPLICATES }       from '../../modules/nf-core/software/gatk4/markduplicates/main'      addParams(options: params.markduplicates_options)
include { GATK4_MARKDUPLICATES_SPARK } from '../../modules/nf-core/software/gatk4/markduplicatesspark/main' addParams(options: params.markduplicates_options)

workflow MARKDUPLICATES {
    take:
        bam_mapped     // channel: [mandatory] meta, bam, bai
        use_gatk_spark //   value: [mandatory] use gatk spark
        save_metrics   //   value: [mandatory] save metrics

    main:

    report_markduplicates = Channel.empty()

    if (use_gatk_spark) {
        GATK4_MARKDUPLICATES_SPARK(bam_mapped, save_metrics)
        report_markduplicates = GATK4_MARKDUPLICATES_SPARK.out.metrics
        bam_markduplicates    = GATK4_MARKDUPLICATES_SPARK.out.bam
    } else {
        GATK4_MARKDUPLICATES(bam_mapped, save_metrics)
        report_markduplicates = GATK4_MARKDUPLICATES.out.metrics
        bam_markduplicates    = GATK4_MARKDUPLICATES.out.bam
    }

    emit:
        bam    = bam_markduplicates
        report = report_markduplicates
}
