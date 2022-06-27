//
// RECALIBRATE_CSV
//

workflow RECALIBRATE_CSV {
    take:
        cram_recalibrated_index // channel: [mandatory] meta, cram, crai

    main:
        // Creating csv files to restart from this step
        cram_recalibrated_index.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/csv") { meta, file, index ->
            patient = meta.patient
            sample  = meta.sample
            gender  = meta.gender
            status  = meta.status
            file = "${params.outdir}/preprocessing/${sample}/recalibrated/${file.name}"
            index = "${params.outdir}/preprocessing/${sample}/recalibrated/${index.name}"
            ["recalibrated.csv", "patient,gender,status,sample,cram,crai\n${patient},${gender},${status},${sample},${file},${index}\n"]
        }
}
