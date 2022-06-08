//
// RECALIBRATE_CSV
//

workflow RECALIBRATE_CSV {
    take:
        cram_recalibrated_index // channel: [mandatory] meta, cram, crai

    main:
        // Creating csv files to restart from this step
        cram_recalibrated_index.collectFile(storeDir: "${params.outdir}/preprocessing/csv") { meta, cram, crai ->
            patient = meta.patient
            sample  = meta.sample
            gender  = meta.gender
            status  = meta.status
            cram = "${params.outdir}/preprocessing/${sample}/recalibrated/${cram.name}"
            crai = "${params.outdir}/preprocessing/${sample}/recalibrated/${crai.name}"
            ["recalibrated_${sample}.csv", "patient,gender,status,sample,cram,crai\n${patient},${gender},${status},${sample},${cram},${crai}\n"]
        }.collectFile(name: 'recalibrated.csv', keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/preprocessing/csv")
}
