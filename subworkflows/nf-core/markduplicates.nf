/*
========================================================================================
    MARKDUPLICATES
========================================================================================
*/

params.markduplicates_options = [:]

include { GATK4_MARKDUPLICATES }                                from '../../modules/nf-core/software/gatk4/markduplicates/main'             addParams(options: params.markduplicates_options)
include { GATK4_MARKDUPLICATES_SPARK }                          from '../../modules/nf-core/software/gatk4/markduplicatesspark/main'        addParams(options: params.markduplicatesspark_options)
include { GATK4_ESTIMATELIBRARYCOMPLEXITY }                     from '../../modules/nf-core/software/gatk4/estimatelibrarycomplexity/main'  addParams(options: params.estimatelibrarycomplexity_options)
workflow MARKDUPLICATES {
    take:
        bam_mapped     // channel: [mandatory] meta, bam, bai
        use_gatk_spark // value: [mandatory] use gatk spark
        save_metrics   // value: [mandatory] save metrics
        fasta          // channel: [mandatory] fasta
        fai            // channel: [mandatory] fai
        dict           // channel: [mandatory] dict
    main:

    report_markduplicates = Channel.empty()

    if (use_gatk_spark) {
        GATK4_MARKDUPLICATES_SPARK(bam_mapped, fasta, fai, dict)
        GATK4_ESTIMATELIBRARYCOMPLEXITY(bam_mapped, fasta, fai, dict)
        report_markduplicates = GATK4_ESTIMATELIBRARYCOMPLEXITY.out.metrics //Here it will Estimate Library complexity
        bam_markduplicates    = GATK4_MARKDUPLICATES_SPARK.out.cram
    } else {
        GATK4_MARKDUPLICATES(bam_mapped, save_metrics)
        report_markduplicates = GATK4_MARKDUPLICATES.out.metrics
        bam_markduplicates    = GATK4_MARKDUPLICATES.out.bam
    }

    emit:
        bam    = bam_markduplicates
        report = report_markduplicates
}
