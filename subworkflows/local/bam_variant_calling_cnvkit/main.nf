//
// CNVKIT calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { CNVKIT_BATCH } from '../../../modules/nf-core/cnvkit/batch/main'
include { CNVKIT_GENEMETRICS } from '../../../modules/nf-core/cnvkit/genemetrics/main'

workflow BAM_VARIANT_CALLING_CNVKIT {
    take:
    cram                // channel: [mandatory] cram
    fasta               // channel: [mandatory] fasta
    fasta_fai           // channel: [optional]  fasta_fai
    targets             // channel: [mandatory] bed
    reference           // channel: [] cnn

    main:
    generate_pon = false

    CNVKIT_BATCH(cram, fasta, fasta_fai, targets, reference, generate_pon)
    CNVKIT_GENEMETRICS(CNVKIT_BATCH.out.cnr, CNVKIT_BATCH.out.cns)

    versions = CNVKIT_BATCH.out.versions

    emit:
    versions // channel: [ versions.yml ]
}
