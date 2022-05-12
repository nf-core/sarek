//
// PREPARE_RECALIBRATION_CSV
//

workflow PREPARE_RECALIBRATION_CSV {
    take:
        cram_table_bqsr // channel: [mandatory] meta, cram, crai, table

    main:
        // Creating csv files to restart from this step
        cram_table_bqsr.collectFile(storeDir: "${params.outdir}/preprocessing/csv") { meta, cram, crai, table ->
            patient = meta.patient
            sample  = meta.sample
            gender  = meta.gender
            status  = meta.status
            cram = "${params.outdir}/preprocessing/${sample}/markduplicates/${sample}.md.cram"
            crai = "${params.outdir}/preprocessing/${sample}/markduplicates/${sample}.md.cram.crai"
            table = "${params.outdir}/preprocessing/${sample}/recal_table/${sample}.recal.table"
            ["markduplicates_${sample}.csv", "patient,gender,status,sample,cram,crai,table\n${patient},${gender},${status},${sample},${cram},${crai},${table}\n"]
        }.collectFile(name: 'markduplicates.csv', keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/preprocessing/csv")
}
