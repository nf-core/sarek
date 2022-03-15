//
// MAPPING_CSV
//

workflow MAPPING_CSV {
    take:
        bam_indexed         // channel: [mandatory] meta, bam, bai

    main:
        csv_bam_mapped = bam_indexed.map { meta, bam, bai -> [meta] }
        // Creating csv files to restart from this step
        csv_bam_mapped.collectFile(storeDir: "${params.outdir}/preprocessing/csv") { meta ->
            patient = meta.patient[0]
            sample  = meta.sample[0]
            gender  = meta.gender[0]
            status  = meta.status[0]
            bam   = "${params.outdir}/preprocessing/${sample}/mapped/${sample}.bam"
            bai   = "${params.outdir}/preprocessing/${sample}/mapped/${sample}.bam.bai"
            ["mapped_${sample}.csv", "patient,gender,status,sample,bam,bai\n${patient},${gender},${status},${sample},${bam},${bai}\n"]
        }.collectFile(name: "mapped.csv", keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/preprocessing/csv")
}
