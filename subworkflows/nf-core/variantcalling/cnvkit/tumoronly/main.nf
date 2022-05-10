//
// CNV calling TUMOR_ONLY
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include {CNVKIT_BATCH as CNVKIT_BATCH_TUMORONLY   } from '../../../../../modules/nf-core/modules/cnvkit/batch/main'

workflow RUN_CNVKIT_TUMORONLY {
    take:
        cram_recalibrated   // channel: [mandatory] cram tumor
        fasta               // channel: [mandatory] fasta
        targets             // channel: [mandatory] bed
        reference           // channel: [mandatory] cnn

    main:
        ch_versions = Channel.empty()

        CNVKIT_BATCH_TUMORONLY(cram_recalibrated, fasta, targets, reference)

        ch_versions = ch_versions.mix(CNVKIT_BATCH_TUMORONLY.out.versions.first())

    emit:
        versions = ch_versions         // channel: [ versions.yml ]

}
