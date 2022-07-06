//
// PREPARE_RECALIBRATION_CSV
//

workflow PREPARE_RECALIBRATION_CSV {
    take:
        cram_table_bqsr // channel: [mandatory] meta, cram, crai, table

    main:
        // Creating csv files to restart from this step
        cram_table_bqsr.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/csv") { meta, cram, crai, table ->
            patient = meta.patient
            sample  = meta.sample
            sex     = meta.sex
            status  = meta.status
            suffix_aligned = params.save_output_as_bam ? "bam" : "cram"
            suffix_index   = params.save_output_as_bam ? "bam.bai" : "cram.crai"
            cram = "${params.outdir}/preprocessing/${sample}/markduplicates/${cram.baseName}.${suffix_aligned}"
            crai = "${params.outdir}/preprocessing/${sample}/markduplicates/${crai.baseName.minus(".cram")}.${suffix_index}"
            table = "${params.outdir}/preprocessing/${sample}/recal_table/${sample}.recal.table"
            ["markduplicates.csv", "patient,sex,status,sample,cram,crai,table\n${patient},${sex},${status},${sample},${cram},${crai},${table}\n"]
        }
}
