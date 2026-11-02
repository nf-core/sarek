//
// CHANNEL_VARIANT_CALLING_CREATE_CSV
//

workflow CHANNEL_VARIANT_CALLING_CREATE_CSV {
    take:
        vcf_to_annotate         // channel: [mandatory] meta, vcf, publish_dir
        outdir                  //

    main:
        // Creating csv files to restart from this step
        vcf_to_annotate.collectFile(keepHeader: true, skip: 1,sort: true, storeDir: "${outdir}/csv"){ meta, vcf, publish_dir ->
            def patient       = meta.patient
            def sample        = meta.id
            def variantcaller = meta.variantcaller
            vcf = "${outdir}/variant_calling/${publish_dir}/${meta.id}/${vcf.getName()}"
            ["variantcalled.csv", "patient,sample,variantcaller,vcf\n${patient},${sample},${variantcaller},${vcf}\n"]
        }
}
