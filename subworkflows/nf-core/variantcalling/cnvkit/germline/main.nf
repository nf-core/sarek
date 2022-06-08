//
// CNV calling GERMLINE
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { CNVKIT_BATCH as CNVKIT_BATCH_GERMLINE                         } from '../../../../../modules/nf-core/modules/cnvkit/batch/main'

workflow RUN_CNVKIT_GERMLINE {
    take:
        cram_recalibrated   // channel: [mandatory] cram
        fasta               // channel: [mandatory] fasta
        fasta_fai           // channel: [optional]  fasta_fai
        targets             // channel: [mandatory] bed
        reference           // channel: [] cnn

    main:
        ch_versions = Channel.empty()

        CNVKIT_BATCH_GERMLINE(cram_recalibrated, fasta, fasta_fai, targets, [])

        ch_versions = ch_versions.mix(CNVKIT_BATCH_GERMLINE.out.versions)

    emit:
        versions = ch_versions              // channel: [ versions.yml ]

}
