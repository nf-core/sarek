/*
========================================================================================
    RECALIBRATE_CSV
========================================================================================
*/

workflow RECALIBRATE_CSV {
    take:
        bam_recalibrated_index // channel: [mandatory] meta, bam, bai

    main:
        // Creating csv files to restart from this step
        bam_recalibrated_index.collectFile(storeDir: "${params.outdir}/preprocessing/csv") { meta, bam, bai ->
            patient = meta.patient
            sample  = meta.sample
            gender  = meta.gender
            status  = meta.status
            bam = "${params.outdir}/preprocessing/${sample}/recalibrated/${sample}.recal.bam"
            bai = "${params.outdir}/preprocessing/${sample}/recalibrated/${sample}.recal.bam.bai"
            ["recalibrated_${sample}.csv", "patient,gender,status,sample,bam,bai\n${patient},${gender},${status},${sample},${bam},${bai}\n"]
        }.collectFile(name: 'recalibrated.csv', keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/preprocessing/csv")
}
