/*
========================================================================================
    MARKDUPLICATES_CSV
========================================================================================
*/

workflow MARKDUPLICATES_CSV {
    take:
        bam_markduplicates // channel: [mandatory] meta, bam, bai

    main:
        // Creating csv files to restart from this step
        bam_markduplicates.collectFile(storeDir: "${params.outdir}/preprocessing/csv") { meta, bam, bai ->
            patient = meta.patient
            sample  = meta.sample
            gender  = meta.gender
            status  = meta.status
            bam   = "${params.outdir}/preprocessing/${sample}/markduplicates/${sample}.md.bam"
            bai   = "${params.outdir}/preprocessing/${sample}/markduplicates/${sample}.md.bam.bai"
            table = "${params.outdir}/preprocessing/${sample}/markduplicates/${sample}.recal.table"
            ["markduplicates_${sample}.csv", "patient,gender,status,sample,bam,bai,table\n${patient},${gender},${status},${sample},${bam},${bai},${table}\n"]
        }.collectFile(name: 'markduplicates.csv', keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/preprocessing/csv")
}
