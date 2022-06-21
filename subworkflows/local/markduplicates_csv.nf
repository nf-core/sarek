//
// MARKDUPLICATES_CSV
//

workflow MARKDUPLICATES_CSV {
    take:
        cram_markduplicates // channel: [mandatory] meta, cram, crai

    main:
        // Creating csv files to restart from this step
        cram_markduplicates.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/csv") { meta, file, index ->
            patient        = meta.patient
            sample         = meta.sample
            gender         = meta.gender
            status         = meta.status
            suffix_aligned = params.save_output_as_bam ? "bam" : "cram"
            suffix_index   = params.save_output_as_bam ? "bam.bai" : "cram.crai"
            file   = "${params.outdir}/preprocessing/${sample}/markduplicates/${file.baseName}.${suffix_aligned}"
            index   = "${params.outdir}/preprocessing/${sample}/markduplicates/${index.baseName.minus(".cram")}.${suffix_index}"
            ["markduplicates_no_table.csv", "patient,gender,status,sample,cram,crai\n${patient},${gender},${status},${sample},${file},${index}\n"]
        }
}
