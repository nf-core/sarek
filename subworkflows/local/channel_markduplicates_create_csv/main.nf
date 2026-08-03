//
// CHANNEL_MARKDUPLICATES_CREATE_CSV
//

workflow CHANNEL_MARKDUPLICATES_CREATE_CSV {
    take:
        cram_markduplicates     // channel: [mandatory] meta, file, index
        csv_subfolder           //
        outdir                  //

    main:
        // Creating csv files to restart from this step
        cram_markduplicates.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${outdir}/csv") { meta, file, index ->
            def patient      = meta.patient
            def sample       = meta.sample
            def sex          = meta.sex
            def status       = meta.status
            def is_bam       = file.name.endsWith('.bam')
            def type         = is_bam ? "bam" : "cram"
            def type_index   = is_bam ? "bai" : "crai"
            // In BAM mode the upstream index is Picard's *.md.bai (from --CREATE_INDEX);
            // GATK4_MARKDUPLICATES publishDir renames it to *.md.bam.bai on disk, so use
            // the standard name in the restart CSV.
            def index_name   = is_bam && index.name.endsWith('.bai') && !index.name.endsWith('.bam.bai')
                ? "${file.name}.bai"
                : index.name
            def align_file   = "${outdir}/preprocessing/${csv_subfolder}/${sample}/${file.name}"
            def align_index  = "${outdir}/preprocessing/${csv_subfolder}/${sample}/${index_name}"

            ["markduplicates_no_table.csv", "patient,sex,status,sample,${type},${type_index}\n${patient},${sex},${status},${sample},${align_file},${align_index}\n"]
        }
}
