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
        cram_recalibrated_index.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${outdir}/csv") { meta, file, index ->
            def patient    = meta.patient
            def sample     = meta.sample
            def sex        = meta.sex
            def status     = meta.status
            def file_path  = "${outdir}/preprocessing/recalibrated/${sample}/${file.name}"
            def index_path = "${outdir}/preprocessing/recalibrated/${sample}/${index.name}"

            def type = save_output_as_bam ? "bam" : "cram"
            def type_index = save_output_as_bam ? "bai" : "crai"

            ["recalibrated.csv", "patient,sex,status,sample,${type},${type_index}\n${patient},${sex},${status},${sample},${file_path},${index_path}\n"]
        }
}
