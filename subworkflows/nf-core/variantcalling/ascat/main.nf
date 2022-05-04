include { ASCAT } from '../../modules/nf-core/modules/ascat/main'

workflow RUN_ASCAT {

    take:
    cram                     // channel: [mandatory] [meta, normal_cram, normal_crai, interval]
    allele_files             // channel: [mandatory]
    loci_files               // channel: [mandatory]

    main:

    ch_versions = Channel.empty()

    ASCAT(cram,allele_files,loci_files)

    emit:
    versions = ch_versions
}
