//
// MAPPING_CSV
//

workflow MAPPING_CSV {
    take:
        bam_indexed         // channel: [mandatory] meta, bam, bai

    main:
        // Creating csv files to restart from this step
        bam_indexed.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/csv") { meta, bam, bai ->
            patient = meta.patient
            sample  = meta.sample
            gender  = meta.gender
            status  = meta.status
            bam   = "${params.outdir}/preprocessing/${sample}/mapped/${bam.name}"
            bai   = "${params.outdir}/preprocessing/${sample}/mapped/${bai.name}"
            ["mapped.csv", "patient,gender,status,sample,bam,bai\n${patient},${gender},${status},${sample},${bam},${bai}\n"]
        }
}
