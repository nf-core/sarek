//
// VARIANTCALLING_CSV
//

workflow VARIANTCALLING_CSV {
    take:
        vcf_to_annotate// channel: [mandatory] meta, vcf

    main:
        // Creating csv files to restart from this step
        //[patient:test3, sample:sample3, status:0, gender:XX, id:sample3, num_intervals:2, variantcaller:Strelka]
        //[[patient:test3, normal_id:sample3, tumor_id:sample4, gender:XX, id:sample4_vs_sample3, num_intervals:2, variantcaller:Strelka, type:indels], /Users/monarchy/Projects/Coding/sarek/work/47/55f3579340cb86b64005cfad606ee0/sample4_vs_sample3.somatic_indels.vcf.gz]
        blub = vcf_to_annotate.collectFile(storeDir: "${params.outdir}/variant_calling/csv") { meta, vcf ->
            patient       = meta.patient
            sample        = meta.id
            gender        = meta.gender
            variantcaller = meta.variantcaller
            vcf = "${params.outdir}/variant_calling/${meta.id}/${variantcaller}/${meta.type}/${vcf.getName()}"
            ["${variantcaller}_${meta.type}_${meta.id}.csv", "patient,gender,sample,variantcaller,type,vcf\n${patient},${gender},${sample},${variantcaller},${meta.type},${vcf}\n"]
        }

        blub.view()

        blub.collectFile(name: 'variantcalled.csv', keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/variant_calling/csv")
}
