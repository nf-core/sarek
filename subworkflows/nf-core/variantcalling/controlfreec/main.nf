include { CONTROLFREEC }    from '../../../../modules/nf-core/modules/controlfreec/main'
include { SAMTOOLS_MPILEUP as MPILEUP }    from '../../../../modules/nf-core/modules/samtools/mpileup/main'

workflow RUN_CONTROLFREEC {
    take:
    cram                     // channel: [mandatory] [meta, cram, crai, interval]
    fasta                    // channel: [mandatory]
    fasta_fai                // channel: [mandatory]
    intervals_bed_gz         // channel: [optional]  Contains a bed.gz file of all intervals combined provided with the cram input(s). Mandatory if interval files are used.
    num_intervals            //     val: [optional]  Number of used intervals, mandatory when intervals are provided.

    main:

    ch_versions = Channel.empty()

    MPILEUP(cram, fasta)
    //MergeMpileup
    //CONTROLFREEC()
    //CONTROLFREECVis


    emit:
    versions = ch_versions
}
