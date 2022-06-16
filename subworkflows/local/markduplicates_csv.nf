//
// MARKDUPLICATES_CSV
//

workflow MARKDUPLICATES_CSV {
    take:
        cram_markduplicates // channel: [mandatory] meta, cram, crai

    main:
        // Creating csv files to restart from this step
        cram_markduplicates.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/csv") { meta, cram, crai ->
            patient = meta.patient
            sample  = meta.sample
            gender  = meta.gender
            status  = meta.status
            cram   = "${params.outdir}/preprocessing/${sample}/markduplicates/${cram.name}"
            crai   = "${params.outdir}/preprocessing/${sample}/markduplicates/${crai.name}"
            ["markduplicates_no_table.csv", "patient,gender,status,sample,cram,crai\n${patient},${gender},${status},${sample},${cram},${crai}\n"]
        }
}
