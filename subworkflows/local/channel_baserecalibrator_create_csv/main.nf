//
// CHANNEL_BASERECALIBRATOR_CREATE_CSV
//

workflow CHANNEL_BASERECALIBRATOR_CREATE_CSV {
    take:
        cram_table_bqsr         // channel: [mandatory] meta, file, index, table
        tools                   //
        skip_tools              //
        outdir                  //

    main:
        // Creating csv files to restart from this step
        if ( tools && tools.split(',').contains('sentieon_dedup') ) {
            cram_table_bqsr.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${outdir}/csv") { meta, file, index, _table ->
                def patient    = meta.patient
                def sample     = meta.sample
                def sex        = meta.sex
                def status     = meta.status
                def is_bam     = file.name.endsWith('.bam')
                def type       = is_bam ? "bam" : "cram"
                def type_index = is_bam ? "bai" : "crai"
                def align_file  = "${outdir}/preprocessing/sentieon_dedup/${sample}/${file.name}"
                def align_index = "${outdir}/preprocessing/sentieon_dedup/${sample}/${index.name}"
                def table_file  = "${outdir}/preprocessing/recal_table/${sample}/${sample}.recal.table"

                ["markduplicates.csv", "patient,sex,status,sample,${type},${type_index},table\n${patient},${sex},${status},${sample},${align_file},${align_index},${table_file}\n"]
            }
        } else if (!(skip_tools && (skip_tools.split(',').contains('markduplicates')))) {
            cram_table_bqsr.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${outdir}/csv") { meta, file, index, _table ->
                def patient    = meta.patient
                def sample     = meta.sample
                def sex        = meta.sex
                def status     = meta.status
                def is_bam     = file.name.endsWith('.bam')
                def type       = is_bam ? "bam" : "cram"
                def type_index = is_bam ? "bai" : "crai"
                def align_file  = "${outdir}/preprocessing/markduplicates/${sample}/${file.name}"
                def align_index = "${outdir}/preprocessing/markduplicates/${sample}/${index.name}"
                def table_file  = "${outdir}/preprocessing/recal_table/${sample}/${sample}.recal.table"

                ["markduplicates.csv", "patient,sex,status,sample,${type},${type_index},table\n${patient},${sex},${status},${sample},${align_file},${align_index},${table_file}\n"]
            }
        } else {
            cram_table_bqsr.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${outdir}/csv") { meta, file, index, _table ->
                def patient    = meta.patient
                def sample     = meta.sample
                def sex        = meta.sex
                def status     = meta.status
                def is_bam     = file.name.endsWith('.bam')
                def type       = is_bam ? "bam" : "cram"
                def type_index = is_bam ? "bai" : "crai"
                def align_file  = "${outdir}/preprocessing/${sample}/mapped/${file.name}"
                def align_index = "${outdir}/preprocessing/${sample}/mapped/${index.name}"
                def table_file  = "${outdir}/preprocessing/${sample}/recal_table/${sample}.recal.table"

                ["sorted.csv", "patient,sex,status,sample,${type},${type_index},table\n${patient},${sex},${status},${sample},${align_file},${align_index},${table_file}\n"]
            }
        }
}
