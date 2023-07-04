//
// CHANNEL_ALIGN_CREATE_CSV
//

workflow CHANNEL_ALIGN_CREATE_CSV {
    take:
        bam_indexed         // channel: [mandatory] meta, bam, bai

    main:
        // Creating csv files to restart from this step
        bam_indexed.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/csv") { meta, bam, bai ->
            patient = meta.patient
            sample  = meta.sample
            sex     = meta.sex
            status  = meta.status
            bam   = "${params.outdir}/preprocessing/mapped/${sample}/${bam.name}"
            bai   = "${params.outdir}/preprocessing/mapped/${sample}/${bai.name}"

            type = params.save_output_as_bam ? "bam" : "cram"
            type_index = params.save_output_as_bam ? "bai" : "crai"

            ["mapped.csv", "patient,sex,status,sample,${type},${type_index}\n${patient},${sex},${status},${sample},${bam},${bai}\n"]
        }
}
