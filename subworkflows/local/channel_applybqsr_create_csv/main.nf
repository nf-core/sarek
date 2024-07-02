//
// CHANNEL_APPLYBQSR_CREATE_CSV
//

workflow CHANNEL_APPLYBQSR_CREATE_CSV {
    take:
        cram_recalibrated_index // channel: [mandatory] meta, cram, crai
        outdir                  //
        save_output_as_bam      //

    main:
        // Creating csv files to restart from this step
        cram_recalibrated_index.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/csv") { meta, file, index ->
            patient = meta.patient
            sample  = meta.sample
            sex     = meta.sex
            status  = meta.status
            file = "${outdir}/preprocessing/recalibrated/${sample}/${file.name}"
            index = "${outdir}/preprocessing/recalibrated/${sample}/${index.name}"

            type = save_output_as_bam ? "bam" : "cram"
            type_index = save_output_as_bam ? "bai" : "crai"

            ["recalibrated.csv", "patient,sex,status,sample,${type},${type_index}\n${patient},${sex},${status},${sample},${file},${index}\n"]
        }
}
