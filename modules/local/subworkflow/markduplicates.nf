/*
================================================================================
                                MAPPING
================================================================================
*/

include { GATK_MARKDUPLICATES } from '../../nf-core/software/gatk/markduplicates'

workflow MARKDUPLICATES {
    take:
        bam_mapped

    main:

    bam_markduplicates    = bam_mapped
    report_markduplicates = Channel.empty()

    if (!params.skip_markduplicates) {
         GATK_MARKDUPLICATES(bam_mapped)
         report_markduplicates = GATK_MARKDUPLICATES.out.report
         bam_markduplicates    = GATK_MARKDUPLICATES.out.bam
         tsv_markduplicates    = GATK_MARKDUPLICATES.out.tsv

        // Creating TSV files to restart from this step
        tsv_markduplicates.collectFile(storeDir: "${params.outdir}/Preprocessing/TSV") { meta ->
            patient = meta.patient
            sample  = meta.sample
            gender  = meta.gender
            status  = meta.status
            bam   = "${params.outdir}/Preprocessing/${sample}/DuplicatesMarked/${sample}.md.bam"
            bai   = "${params.outdir}/Preprocessing/${sample}/DuplicatesMarked/${sample}.md.bam.bai"
            table = "${params.outdir}/Preprocessing/${sample}/DuplicatesMarked/${sample}.recal.table"
            ["duplicates_marked_${sample}.tsv", "${patient}\t${gender}\t${status}\t${sample}\t${bam}\t${bai}\t${table}\n"]
        }

        tsv_markduplicates.map { meta ->
            patient = meta.patient
            sample  = meta.sample
            gender  = meta.gender
            status  = meta.status
            bam   = "${params.outdir}/Preprocessing/${sample}/DuplicatesMarked/${sample}.md.bam"
            bai   = "${params.outdir}/Preprocessing/${sample}/DuplicatesMarked/${sample}.md.bam.bai"
            table = "${params.outdir}/Preprocessing/${sample}/DuplicatesMarked/${sample}.recal.table"
            "${patient}\t${gender}\t${status}\t${sample}\t${bam}\t${bai}\t${table}\n"
        }.collectFile(name: 'duplicates_marked.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV")
    } else {
        tsv_no_markduplicates = bam_markduplicates.map { meta, bam, bai -> [meta] }

        // Creating TSV files to restart from this step
        tsv_no_markduplicates.collectFile(storeDir: "${params.outdir}/Preprocessing/TSV") { meta ->
            patient = meta.patient[0]
            sample  = meta.sample[0]
            gender  = meta.gender[0]
            status  = meta.status[0]
            bam   = "${params.outdir}/Preprocessing/${sample}/Mapped/${sample}.bam"
            bai   = "${params.outdir}/Preprocessing/${sample}/Mapped/${sample}.bam.bai"
            table = "${params.outdir}/Preprocessing/${sample}/Mapped/${sample}.recal.table"
            ["mapped_no_duplicates_marked_${sample}.tsv", "${patient}\t${gender}\t${status}\t${sample}\t${bam}\t${bai}\t${table}\n"]
        }

        tsv_no_markduplicates.map { meta ->
            patient = meta.patient[0]
            sample  = meta.sample[0]
            gender  = meta.gender[0]
            status  = meta.status[0]
            bam   = "${params.outdir}/Preprocessing/${sample}/Mapped/${sample}.bam"
            bai   = "${params.outdir}/Preprocessing/${sample}/Mapped/${sample}.bam.bai"
            table = "${params.outdir}/Preprocessing/${sample}/Mapped/${sample}.recal.table"
            "${patient}\t${gender}\t${status}\t${sample}\t${bam}\t${bai}\t${table}\n"
        }.collectFile(name: 'mapped_no_duplicates_marked.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV")
    }

    emit:
        bam    = bam_markduplicates
        report = report_markduplicates
}
