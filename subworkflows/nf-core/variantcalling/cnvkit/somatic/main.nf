//
// CNV calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { CNVKIT_BATCH     } from '../../../../../modules/nf-core/modules/cnvkit/batch/main'

workflow RUN_CNVKIT_SOMATIC {
    take:
        cram_pair       // channel: [mandatory] cram tumor, cram normal
        fasta           // channel: [mandatory] fasta
        targets         // channel: [mandatory] bed
        reference       // channel: [mandatory]

    main:
        ch_versions = Channel.empty()

        CNVKIT_BATCH(cram_pair, fasta, targets, reference)

        ch_versions = ch_versions.mix(CNVKIT_BATCH.out.versions.first())

    emit:
        versions = ch_versions              // channel: [ versions.yml ]


}
