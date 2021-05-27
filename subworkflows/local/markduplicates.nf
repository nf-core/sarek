/*
================================================================================
                                 MARKDUPLICATES
================================================================================
*/

params.markduplicates_options = [:]

include { GATK4_MARKDUPLICATES }       from '../../modules/nf-core/software/gatk4/markduplicates/main' addParams(options: params.markduplicates_options)
include { GATK4_MARKDUPLICATES_SPARK } from '../../modules/nf-core/software/gatk4/markduplicates/main' addParams(options: params.markduplicates_options)

workflow MARKDUPLICATES {
    take:
        bam_mapped   // channel: [mandatory] bam_mapped
        save_metrics //   value: [mandatory] save metrics

    main:

    bam_markduplicates    = bam_mapped
    report_markduplicates = Channel.empty()

    if (params.use_gatk_spark) {
        GATK4_MARKDUPLICATES_SPARK(bam_mapped, save_metrics)
        report_markduplicates = GATK4_MARKDUPLICATES_SPARK.out.metrics
        bam_markduplicates    = GATK4_MARKDUPLICATES_SPARK.out.bam
    } else {
        GATK4_MARKDUPLICATES(bam_mapped, save_metrics)
        report_markduplicates = GATK4_MARKDUPLICATES.out.metrics
        bam_markduplicates    = GATK4_MARKDUPLICATES.out.bam
    }

    // Creating TSV files to restart from this step
    bam_markduplicates.collectFile(storeDir: "${params.outdir}/preprocessing/tsv") { meta, bam, bai ->
        patient = meta.patient
        sample  = meta.sample
        gender  = meta.gender
        status  = meta.status
        bam   = "${params.outdir}/preprocessing/${sample}/markduplicates/${sample}.md.bam"
        bai   = "${params.outdir}/preprocessing/${sample}/markduplicates/${sample}.md.bam.bai"
        table = "${params.outdir}/preprocessing/${sample}/markduplicates/${sample}.recal.table"
        ["markduplicates_${sample}.tsv", "${patient}\t${gender}\t${status}\t${sample}\t${bam}\t${bai}\t${table}\n"]
    }

    bam_markduplicates.map { meta, bam, bai ->
        patient = meta.patient
        sample  = meta.sample
        gender  = meta.gender
        status  = meta.status
        bam   = "${params.outdir}/preprocessing/${sample}/markduplicates/${sample}.md.bam"
        bai   = "${params.outdir}/preprocessing/${sample}/markduplicates/${sample}.md.bam.bai"
        table = "${params.outdir}/preprocessing/${sample}/markduplicates/${sample}.recal.table"
        "${patient}\t${gender}\t${status}\t${sample}\t${bam}\t${bai}\t${table}\n"
    }.collectFile(name: 'markduplicates.tsv', sort: true, storeDir: "${params.outdir}/preprocessing/tsv")

    emit:
        bam    = bam_markduplicates
        report = report_markduplicates
}
