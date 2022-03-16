//
// PREPARE_RECALIBRATION_CSV
//

workflow PREPARE_RECALIBRATION_CSV {
    take:
        table_bqsr // channel: [mandatory] meta, table

    main:
        // Creating csv files to restart from this step
        table_bqsr.collectFile(storeDir: "${params.outdir}/preprocessing/csv") { meta, table ->
            patient = meta.patient
            sample  = meta.sample
            gender  = meta.gender
            status  = meta.status
            cram = "${params.outdir}/preprocessing/${sample}/markduplicates/${sample}.md.cram"
            crai = "${params.outdir}/preprocessing/${sample}/markduplicates/${sample}.md.cram.crai"
            ["markduplicates_no_table_${sample}.csv", "patient,gender,status,sample,cram,crai\n${patient},${gender},${status},${sample},${cram},${crai}\n"]
        }.collectFile(name: 'markduplicates_no_table.csv', keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/preprocessing/csv")
}
