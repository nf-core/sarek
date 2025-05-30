//
// CHANNEL_BASERECALIBRATOR_CREATE_CSV
//

workflow CHANNEL_BASERECALIBRATOR_CREATE_CSV {
    take:
        cram_table_bqsr         // channel: [mandatory] meta, cram, crai, table
        tools                   //
        skip_tools              //
        outdir                  //
        save_output_as_bam      //

    main:
        // Creating csv files to restart from this step
        if ( tools && tools.split(',').contains('sentieon_dedup') ) {
            cram_table_bqsr.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${outdir}/csv") { meta, cram, crai, table ->

                def patient        = meta.patient
                def sample         = meta.sample
                def sex            = meta.sex
                def status         = meta.status
                def suffix_aligned = save_output_as_bam ? "bam" : "cram"
                def suffix_index   = save_output_as_bam ? "bam.bai" : "cram.crai"
                def file_path      = "${outdir}/preprocessing/sentieon_dedup/${sample}/${cram.baseName}.${suffix_aligned}"
                def index_path     = "${outdir}/preprocessing/sentieon_dedup/${sample}/${crai.baseName.minus(".cram")}.${suffix_index}"
                def table_path     = "${outdir}/preprocessing/recal_table/${sample}/${sample}.recal.table"

                def type = save_output_as_bam ? "bam" : "cram"
                def type_index = save_output_as_bam ? "bai" : "crai"

                ["markduplicates.csv", "patient,sex,status,sample,${type},${type_index},table\n${patient},${sex},${status},${sample},${file_path},${index_path},${table_path}\n"]
            }
        } else if (!(skip_tools && (skip_tools.split(',').contains('markduplicates')))) {
            cram_table_bqsr.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${outdir}/csv") { meta, cram, crai, table ->

                def patient        = meta.patient
                def sample         = meta.sample
                def sex            = meta.sex
                def status         = meta.status
                def suffix_aligned = save_output_as_bam ? "bam" : "cram"
                def suffix_index   = save_output_as_bam ? "bam.bai" : "cram.crai"
                def file_path      = "${outdir}/preprocessing/markduplicates/${sample}/${cram.baseName}.${suffix_aligned}"
                def index_path     = "${outdir}/preprocessing/markduplicates/${sample}/${crai.baseName.minus(".cram")}.${suffix_index}"
                def table_path     = "${outdir}/preprocessing/recal_table/${sample}/${sample}.recal.table"

                def type = save_output_as_bam ? "bam" : "cram"
                def type_index = save_output_as_bam ? "bai" : "crai"

                ["markduplicates.csv", "patient,sex,status,sample,${type},${type_index},table\n${patient},${sex},${status},${sample},${file_path},${index_path},${table_path}\n"]
            }
        } else {
            cram_table_bqsr.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${outdir}/csv") { meta, cram, crai, table ->
                def patient        = meta.patient
                def sample         = meta.sample
                def sex            = meta.sex
                def status         = meta.status
                def suffix_aligned = save_output_as_bam ? "bam" : "cram"
                def suffix_index   = save_output_as_bam ? "bam.bai" : "cram.crai"
                def file_path      = "${outdir}/preprocessing/${sample}/mapped/${cram.baseName}.${suffix_aligned}"
                def index_path     = "${outdir}/preprocessing/${sample}/mapped/${crai.baseName.minus(".cram")}.${suffix_index}"
                def table_path     = "${outdir}/preprocessing/${sample}/recal_table/${sample}.recal.table"

                def type = save_output_as_bam ? "bam" : "cram"
                def type_index = save_output_as_bam ? "bai" : "crai"

                ["sorted.csv", "patient,sex,status,sample,${type},${type_index},table\n${patient},${sex},${status},${sample},${file_path},${index_path},${table_path}\n"]
            }
        }
}
