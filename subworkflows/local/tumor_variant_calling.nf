//
// TUMOR VARIANT CALLING
//



workflow TUMOR_VARIANT_CALLING {
    take:
        tools
        cram              // channel: [mandatory] cram
        dbsnp             // channel: [mandatory] dbsnp
        dbsnp_tbi         // channel: [mandatory] dbsnp_tbi
        dict              // channel: [mandatory] dict
        fasta             // channel: [mandatory] fasta
        fasta_fai         // channel: [mandatory] fasta_fai
        intervals         // channel: [mandatory] intervals
        num_intervals
        target_bed        // channel: [optional]  target_bed
        target_bed_gz_tbi // channel: [optional]  target_bed_gz_tbi

    main:


    emit:

}
