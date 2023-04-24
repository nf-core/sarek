//
// CHANNEL_MARKDUPLICATES_CREATE_CSV
//

workflow CHANNEL_MARKDUPLICATES_CREATE_CSV {
    take:
        cram_markduplicates // channel: [mandatory] meta, cram, crai
        csv_subfolder
        outdir
        save_output_as_bam

    main:
        // Creating csv files to restart from this step
        cram_markduplicates.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${outdir}/csv") { meta, file, index ->
            patient        = meta.patient
            sample         = meta.sample
            sex            = meta.sex
            status         = meta.status
            suffix_aligned = save_output_as_bam ? "bam" : "cram"
            suffix_index   = save_output_as_bam ? "bam.bai" : "cram.crai"
            file   = "${outdir}/preprocessing/${csv_subfolder}/${sample}/${file.baseName}.${suffix_aligned}"
            index   = "${outdir}/preprocessing/${csv_subfolder}/${sample}/${index.baseName.minus(".cram")}.${suffix_index}"

            type = save_output_as_bam ? "bam" : "cram"
            type_index = save_output_as_bam ? "bai" : "crai"

            ["markduplicates_no_table.csv", "patient,sex,status,sample,${type},${type_index}\n${patient},${sex},${status},${sample},${file},${index}\n"]
        }
}
