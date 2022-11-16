//
// CNVKit calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { CNVKIT_BATCH } from '../../../modules/nf-core/cnvkit/batch/main'

workflow BAM_VARIANT_CALLING_CNVKIT {
    take:
        cram_recalibrated   // channel: [mandatory] cram
        fasta               // channel: [mandatory] fasta
        fasta_fai           // channel: [optional]  fasta_fai
        targets             // channel: [mandatory] bed
        reference           // channel: [] cnn

    main:
        ch_versions = Channel.empty()

        CNVKIT_BATCH(cram_recalibrated, fasta, fasta_fai, targets, reference)

        ch_versions = ch_versions.mix(CNVKIT_BATCH.out.versions)

    emit:
        versions = ch_versions              // channel: [ versions.yml ]

}
