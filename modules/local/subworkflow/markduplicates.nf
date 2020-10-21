/*
================================================================================
                                 MARKDUPLICATES
================================================================================
*/

params.adam_transformalignments_options = [:]
params.markduplicates_options           = [:]
params.samtools_index_options           = [:]

include { ADAM_TRANSFORMALIGNMENTS } from '../../nf-core/software/adam/transformalignments' addParams(options: params.adam_transformalignments_options)
include { GATK_MARKDUPLICATES }      from '../../nf-core/software/gatk/markduplicates'      addParams(options: params.markduplicates_options)
include { SAMTOOLS_INDEX }           from '../../nf-core/software/samtools/index'           addParams(options: params.samtools_index_options)

workflow MARKDUPLICATES {
    take:
        bam_mapped // channel: [mandatory] bam_mapped
        step       //   value: [mandatory] starting step

    main:

    bam_markduplicates    = bam_mapped
    report_markduplicates = Channel.empty()

    if (step == "mapping") {
        if (!params.skip_markduplicates) {
            if (params.markduplicates == 'gatk') {
                GATK_MARKDUPLICATES(bam_mapped)

                bam_markduplicates    = GATK_MARKDUPLICATES.out.bam
                report_markduplicates = GATK_MARKDUPLICATES.out.report
                tsv_markduplicates    = GATK_MARKDUPLICATES.out.tsv
            } else if (params.markduplicates == 'adam') {
                ADAM_TRANSFORMALIGNMENTS(bam_mapped)

                bam_markduplicates    = SAMTOOLS_INDEX(ADAM_TRANSFORMALIGNMENTS.out.bam)
                tsv_markduplicates    = ADAM_TRANSFORMALIGNMENTS.out.tsv
            }

            // Creating TSV files to restart from this step
            tsv_markduplicates.collectFile(storeDir: "${params.outdir}/preprocessing/tsv") { meta ->
                patient = meta.patient
                sample  = meta.sample
                gender  = meta.gender
                status  = meta.status
                bam   = "${params.outdir}/preprocessing/${sample}/markduplicates/${sample}.md.bam"
                bai   = "${params.outdir}/preprocessing/${sample}/markduplicates/${sample}.md.bam.bai"
                table = "${params.outdir}/preprocessing/${sample}/markduplicates/${sample}.recal.table"
                ["markduplicates_${sample}.tsv", "${patient}\t${gender}\t${status}\t${sample}\t${bam}\t${bai}\t${table}\n"]
            }

            tsv_markduplicates.map { meta ->
                patient = meta.patient
                sample  = meta.sample
                gender  = meta.gender
                status  = meta.status
                bam   = "${params.outdir}/preprocessing/${sample}/markduplicates/${sample}.md.bam"
                bai   = "${params.outdir}/preprocessing/${sample}/markduplicates/${sample}.md.bam.bai"
                table = "${params.outdir}/preprocessing/${sample}/markduplicates/${sample}.recal.table"
                "${patient}\t${gender}\t${status}\t${sample}\t${bam}\t${bai}\t${table}\n"
            }.collectFile(name: 'markduplicates.tsv', sort: true, storeDir: "${params.outdir}/preprocessing/tsv")
        } else {
            tsv_no_markduplicates = bam_markduplicates.map { meta, bam, bai -> [meta] }

            // Creating TSV files to restart from this step
            tsv_no_markduplicates.collectFile(storeDir: "${params.outdir}/preprocessing/tsv") { meta ->
                patient = meta.patient[0]
                sample  = meta.sample[0]
                gender  = meta.gender[0]
                status  = meta.status[0]
                bam   = "${params.outdir}/preprocessing/${sample}/mapped/${sample}.bam"
                bai   = "${params.outdir}/preprocessing/${sample}/mapped/${sample}.bam.bai"
                table = "${params.outdir}/preprocessing/${sample}/mapped/${sample}.recal.table"
                ["mapped_no_markduplicates_${sample}.tsv", "${patient}\t${gender}\t${status}\t${sample}\t${bam}\t${bai}\t${table}\n"]
            }

            tsv_no_markduplicates.map { meta ->
                patient = meta.patient[0]
                sample  = meta.sample[0]
                gender  = meta.gender[0]
                status  = meta.status[0]
                bam   = "${params.outdir}/preprocessing/${sample}/mapped/${sample}.bam"
                bai   = "${params.outdir}/preprocessing/${sample}/mapped/${sample}.bam.bai"
                table = "${params.outdir}/preprocessing/${sample}/mapped/${sample}.recal.table"
                "${patient}\t${gender}\t${status}\t${sample}\t${bam}\t${bai}\t${table}\n"
            }.collectFile(name: 'mapped_no_markduplicates.tsv', sort: true, storeDir: "${params.outdir}/preprocessing/tsv")
        }
    }

    emit:
        bam    = bam_markduplicates
        report = report_markduplicates
}
