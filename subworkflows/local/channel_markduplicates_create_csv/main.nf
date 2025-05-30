//
// CHANNEL_MARKDUPLICATES_CREATE_CSV
//

workflow CHANNEL_MARKDUPLICATES_CREATE_CSV {
    take:
        cram_markduplicates     // channel: [mandatory] meta, cram, crai
        csv_subfolder           //
        outdir                  //
        save_output_as_bam      //

    main:
        // Creating csv files to restart from this step
        cram_markduplicates.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${outdir}/csv") { meta, file, index ->
            def patient        = meta.patient
            def sample         = meta.sample
            def sex            = meta.sex
            def status         = meta.status
            def suffix_aligned = save_output_as_bam ? "bam" : "cram"
            def suffix_index   = save_output_as_bam ? "bam.bai" : "cram.crai"
            def file_path      = "${outdir}/preprocessing/${csv_subfolder}/${sample}/${file.baseName}.${suffix_aligned}"
            def index_path     = "${outdir}/preprocessing/${csv_subfolder}/${sample}/${index.baseName.minus(".cram")}.${suffix_index}"

            def type = save_output_as_bam ? "bam" : "cram"
            def type_index = save_output_as_bam ? "bai" : "crai"

            ["markduplicates_no_table.csv", "patient,sex,status,sample,${type},${type_index}\n${patient},${sex},${status},${sample},${file_path},${index_path}\n"]
        }
}
