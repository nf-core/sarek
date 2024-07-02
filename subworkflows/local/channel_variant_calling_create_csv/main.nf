//
// CHANNEL_VARIANT_CALLING_CREATE_CSV
//

workflow CHANNEL_VARIANT_CALLING_CREATE_CSV {
    take:
        vcf_to_annotate         // channel: [mandatory] meta, vcf
        outdir                  //

    main:
        // Creating csv files to restart from this step
        vcf_to_annotate.collectFile(keepHeader: true, skip: 1,sort: true, storeDir: "${outdir}/csv"){ meta, vcf ->
            patient       = meta.patient
            sample        = meta.id
            variantcaller = meta.variantcaller
            vcf = "${outdir}/variant_calling/${variantcaller}/${meta.id}/${vcf.getName()}"
            ["variantcalled.csv", "patient,sample,variantcaller,vcf\n${patient},${sample},${variantcaller},${vcf}\n"]
        }
}
