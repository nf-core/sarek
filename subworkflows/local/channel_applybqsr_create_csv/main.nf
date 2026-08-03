//
// CHANNEL_APPLYBQSR_CREATE_CSV
//

workflow CHANNEL_APPLYBQSR_CREATE_CSV {
    take:
    cram_recalibrated_index // channel: [mandatory] meta, file, index
    outdir //

    main:
    // Creating csv files to restart from this step
    cram_recalibrated_index.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${outdir}/csv") { meta, file, index ->
        def patient = meta.patient
        def sample = meta.sample
        def sex = meta.sex
        def status = meta.status
        def is_bam = file.name.endsWith('.bam')
        def type = is_bam ? "bam" : "cram"
        def type_index = is_bam ? "bai" : "crai"
        def out_file = "${outdir}/preprocessing/recalibrated/${sample}/${file.name}"
        def out_index = "${outdir}/preprocessing/recalibrated/${sample}/${index.name}"

        ["recalibrated.csv", "patient,sex,status,sample,${type},${type_index}\n${patient},${sex},${status},${sample},${out_file},${out_index}\n"]
    }
}
