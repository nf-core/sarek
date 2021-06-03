/*
========================================================================================
    MAPPING_CSV
========================================================================================
*/

workflow MAPPING_CSV {
    take:
        bam_mapped          // channel: [mandatory] meta, bam, bai
        save_bam_mapped     // boolean: [mandatory] save_bam_mapped
        skip_markduplicates // boolean: [mandatory] skip_markduplicates

    main:
        if (save_bam_mapped) {
            csv_bam_mapped = bam_mapped.map { meta, bam, bai -> [meta] }
            // Creating csv files to restart from this step
            csv_bam_mapped.collectFile(storeDir: "${params.outdir}/preprocessing/csv") { meta ->
                patient = meta.patient[0]
                sample  = meta.sample[0]
                gender  = meta.gender[0]
                status  = meta.status[0]
                bam   = "${params.outdir}/preprocessing/${sample}/mapped/${sample}.bam"
                bai   = "${params.outdir}/preprocessing/${sample}/mapped/${sample}.bam.bai"
                ["mapped_${sample}.csv", "${patient},${gender},${status},${sample},${bam},${bai}\n"]
            }

            csv_bam_mapped.map { meta ->
                patient = meta.patient[0]
                sample  = meta.sample[0]
                gender  = meta.gender[0]
                status  = meta.status[0]
                bam   = "${params.outdir}/preprocessing/${sample}/mapped/${sample}.bam"
                bai   = "${params.outdir}/preprocessing/${sample}/mapped/${sample}.bam.bai"
                "${patient},${gender},${status},${sample},${bam},${bai}\n"
            }.collectFile(name: "mapped.csv", sort: true, storeDir: "${params.outdir}/preprocessing/csv")
        }

        if (skip_markduplicates) {
            csv_bam_mapped = bam_mapped.map { meta, bam, bai -> [meta] }
            // Creating csv files to restart from this step
            csv_bam_mapped.collectFile(storeDir: "${params.outdir}/preprocessing/csv") { meta ->
                patient = meta.patient[0]
                sample  = meta.sample[0]
                gender  = meta.gender[0]
                status  = meta.status[0]
                bam   = "${params.outdir}/preprocessing/${sample}/mapped/${sample}.bam"
                bai   = "${params.outdir}/preprocessing/${sample}/mapped/${sample}.bam.bai"
                table = "${params.outdir}/preprocessing/${sample}/mapped/${sample}.recal.table"
                ["mapped_no_markduplicates_${sample}.csv", "${patient},${gender},${status},${sample},${bam},${bai},${table}\n"]
            }

            csv_bam_mapped.map { meta ->
                patient = meta.patient[0]
                sample  = meta.sample[0]
                gender  = meta.gender[0]
                status  = meta.status[0]
                bam   = "${params.outdir}/preprocessing/${sample}/mapped/${sample}.bam"
                bai   = "${params.outdir}/preprocessing/${sample}/mapped/${sample}.bam.bai"
                table = "${params.outdir}/preprocessing/${sample}/mapped/${sample}.recal.table"
                "${patient},${gender},${status},${sample},${bam},${bai},${table}\n"
            }.collectFile(name: 'mapped_no_markduplicates.csv', sort: true, storeDir: "${params.outdir}/preprocessing/csv")
        }
}
