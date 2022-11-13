//
// CHANNEL_BASERECALIBRATOR_CREATE_CSV
//

workflow CHANNEL_BASERECALIBRATOR_CREATE_CSV {
    take:
        cram_table_bqsr // channel: [mandatory] meta, cram, crai, table
        skip_tools

    main:
        // Creating csv files to restart from this step
        if (!(skip_tools && (skip_tools.split(',').contains('markduplicates')))) {
            cram_table_bqsr.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/csv") { meta, cram, crai, table ->

                patient = meta.patient
                sample  = meta.sample
                sex     = meta.sex
                status  = meta.status
                suffix_aligned = params.save_output_as_bam ? "bam" : "cram"
                suffix_index   = params.save_output_as_bam ? "bam.bai" : "cram.crai"
                cram = "${params.outdir}/preprocessing/markduplicates/${sample}/${cram.baseName}.${suffix_aligned}"
                crai = "${params.outdir}/preprocessing/markduplicates/${sample}/${crai.baseName.minus(".cram")}.${suffix_index}"
                table = "${params.outdir}/preprocessing/recal_table/${sample}/${sample}.recal.table"

                type = params.save_output_as_bam ? "bam" : "cram"
                type_index = params.save_output_as_bam ? "bai" : "crai"

                ["markduplicates.csv", "patient,sex,status,sample,${type},${type_index},table\n${patient},${sex},${status},${sample},${cram},${crai},${table}\n"]
            }
        } else {
            cram_table_bqsr.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/csv") { meta, cram, crai, table ->
                patient = meta.patient
                sample  = meta.sample
                sex     = meta.sex
                status  = meta.status
                suffix_aligned = params.save_output_as_bam ? "bam" : "cram"
                suffix_index   = params.save_output_as_bam ? "bam.bai" : "cram.crai"
                cram = "${params.outdir}/preprocessing/${sample}/mapped/${cram.baseName}.${suffix_aligned}"
                crai = "${params.outdir}/preprocessing/${sample}/mapped/${crai.baseName.minus(".cram")}.${suffix_index}"
                table = "${params.outdir}/preprocessing/${sample}/recal_table/${sample}.recal.table"

                type = params.save_output_as_bam ? "bam" : "cram"
                type_index = params.save_output_as_bam ? "bai" : "crai"

                ["sorted.csv", "patient,sex,status,sample,${type},${type_index},table\n${patient},${sex},${status},${sample},${cram},${crai},${table}\n"]
            }
        }
}
