//
// VARIANTCALLING_CSV
//

workflow VARIANTCALLING_CSV {
    take:
        vcf_to_annotate // channel: [mandatory] meta, vcf

    main:
        // Creating csv files to restart from this step
        vcf_to_annotate.collectFile(storeDir: "${params.outdir}/variant_calling/csv") { meta, vcf ->
            patient       = meta.patient
            sample        = meta.id
            gender        = meta.gender
            variantcaller = meta.variantcaller
            vcf = "${params.outdir}/variant_calling/${meta.id}/${variantcaller}/${vcf.getName()}"
            ["${variantcaller}_${meta.id}.csv", "patient,gender,sample,variantcaller,vcf\n${patient},${gender},${sample},${variantcaller},${vcf}\n"]
        }.collectFile(name: 'variantcalled.csv', keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/variant_calling/csv")
}
