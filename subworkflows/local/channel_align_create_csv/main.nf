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
        bam_indexed.collectFile(keepHeader: true, skip: 1, sort: true) { meta, bam, bai ->
            patient = meta.patient
            sample  = meta.sample
            sex     = meta.sex
            status  = meta.status
            bam   = "${outdir}/preprocessing/mapped/${sample}/${bam.name}"
            bai   = "${outdir}/preprocessing/mapped/${sample}/${bai.name}"

            type = save_output_as_bam ? "bam" : "cram"
            type_index = save_output_as_bam ? "bai" : "crai"

            ["mapped.csv", "patient,sex,status,sample,${type},${type_index}\n${patient},${sex},${status},${sample},${bam},${bai}\n"]
        }
}
