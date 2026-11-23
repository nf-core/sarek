//
// CHANNEL_BASERECALIBRATOR_CREATE_CSV
//

workflow CHANNEL_BASERECALIBRATOR_CREATE_CSV {
    take:
        cram_table_bqsr         // channel: [mandatory] meta, file, index, table
        tools                   //
        skip_tools              //
        step                    //
        sentieon_consensus      //
        outdir                  //

    main:
        // Creating csv files to restart from this step
        if ( 'sentieon_dedup' in tools ) {
            cram_table_bqsr.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${outdir}/csv") { meta, file, index, _table ->
                def patient    = meta.patient
                def sample     = meta.sample
                def sex        = meta.sex
                def status     = meta.status
                def is_bam     = file.name.endsWith('.bam')
                def type       = is_bam ? "bam" : "cram"
                def type_index = is_bam ? "bai" : "crai"
                def dedup_dir   = sentieon_consensus ? "sentieon_consensus" : "sentieon_dedup"
                def align_file  = "${outdir}/preprocessing/${dedup_dir}/${sample}/${file.name}"
                def align_index = "${outdir}/preprocessing/${dedup_dir}/${sample}/${index.name}"
                def table_file  = "${outdir}/preprocessing/recal_table/${sample}/${sample}.recal.table"

                ["markduplicates.csv", "patient,sex,status,sample,${type},${type_index},table\n${patient},${sex},${status},${sample},${align_file},${align_index},${table_file}\n"]
            }
        } else if ('markduplicates' !in skip_tools) {
            cram_table_bqsr.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${outdir}/csv") { meta, file, index, _table ->
                def patient    = meta.patient
                def sample     = meta.sample
                def sex        = meta.sex
                def status     = meta.status
                def is_bam     = file.name.endsWith('.bam')
                def type       = is_bam ? "bam" : "cram"
                def type_index = is_bam ? "bai" : "crai"
                // Picard's *.md.bai is published as *.md.bam.bai, see CHANNEL_MARKDUPLICATES_CREATE_CSV
                def index_name  = is_bam && index.name.endsWith('.md.bai') ? "${file.name}.bai" : index.name
                // From --step prepare_recalibration, the input files are not published, so list them as given
                def align_file  = step == 'prepare_recalibration' ? file.toUriString() : "${outdir}/preprocessing/markduplicates/${sample}/${file.name}"
                def align_index = step == 'prepare_recalibration' ? index.toUriString() : "${outdir}/preprocessing/markduplicates/${sample}/${index_name}"
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
