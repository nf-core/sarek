//
// CNV calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { CNVKIT_BATCH as CNVKIT_BATCH_SOMATIC } from '../../../../../modules/nf-core/modules/cnvkit/batch/main'

workflow RUN_CNVKIT_SOMATIC {
    take:
        cram_pair       // channel: [mandatory] cram tumor, cram normal
        fasta           // channel: [mandatory] fasta
        fasta_fai       // channel: [optional]  fasta_fai
        targets         // channel: [mandatory] bed
        reference       // channel: [optional]

    main:
        ch_versions = Channel.empty()

        CNVKIT_BATCH_SOMATIC(cram_pair, fasta, fasta_fai, targets, [])

        ch_versions = ch_versions.mix(CNVKIT_BATCH_SOMATIC.out.versions.first())

    emit:
        versions = ch_versions              // channel: [ versions.yml ]


}
