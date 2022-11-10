//
// CHANNEL_APPLYBQSR_CREATE_CSV
//

workflow CHANNEL_APPLYBQSR_CREATE_CSV {
    take:
        cram_recalibrated_index // channel: [mandatory] meta, cram, crai

    main:
        // Creating csv files to restart from this step
        cram_recalibrated_index.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/csv") { meta, file, index ->
            patient = meta.patient
            sample  = meta.sample
            sex     = meta.sex
            status  = meta.status
            file = "${params.outdir}/preprocessing/recalibrated/${sample}/${file.name}"
            index = "${params.outdir}/preprocessing/recalibrated/${sample}/${index.name}"
            ["recalibrated.csv", "patient,sex,status,sample,cram,crai\n${patient},${sex},${status},${sample},${file},${index}\n"]
        }
}
