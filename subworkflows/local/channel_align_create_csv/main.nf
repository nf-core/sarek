//
// CHANNEL_ALIGN_CREATE_CSV
//

workflow CHANNEL_ALIGN_CREATE_CSV {
    take:
        bam_indexed             // channel: [mandatory] meta, bam, bai
        outdir                  //
        save_output_as_bam      //

    main:
        // Creating csv files to restart from this step
        bam_indexed.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${outdir}/csv") { meta, bam, bai ->
            def patient  = meta.patient
            def sample   = meta.sample
            def sex      = meta.sex
            def status   = meta.status
            def bam_path = "${outdir}/preprocessing/mapped/${sample}/${bam.name}"
            def bai_path = "${outdir}/preprocessing/mapped/${sample}/${bai.name}"

            def type = save_output_as_bam ? "bam" : "cram"
            def type_index = save_output_as_bam ? "bai" : "crai"

            ["mapped.csv", "patient,sex,status,sample,${type},${type_index}\n${patient},${sex},${status},${sample},${bam_path},${bai_path}\n"]
        }
}
