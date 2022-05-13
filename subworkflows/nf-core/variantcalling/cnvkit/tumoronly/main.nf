//
// CNV calling TUMOR_ONLY
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run


include {CNVKIT_ANTITARGET                          } from '../../../../../modules/local/cnvkit/antitarget'
include {CNVKIT_REFERENCE                           } from '../../../../../modules/local/cnvkit/reference'
include {CNVKIT_BATCH as CNVKIT_BATCH_TUMORONLY     } from '../../../../../modules/nf-core/modules/cnvkit/batch/main'

workflow RUN_CNVKIT_TUMORONLY {
    take:
        cram_recalibrated   // channel: [mandatory] cram tumor
        fasta               // channel: [mandatory] fasta
        targets             // channel: [mandatory] bed
        reference           // channel: [mandatory] cnn

    main:
        ch_versions = Channel.empty()

        // prepare a reference for tumor_only mode based on target_baits

        CNVKIT_ANTITARGET(targets)

        CNVKIT_REFERENCE(fasta, targets, CNVKIT_ANTITARGET.out.BED)

        // use reference for calling CNVs

        CNVKIT_BATCH_TUMORONLY(cram_recalibrated, [], [], CNVKIT_REFERENCE.out.CNN)

        ch_versions = ch_versions.mix(CNVKIT_ANTITARGET.out.versions)
        ch_versions = ch_versions.mix(CNVKIT_REFERENCE.out.versions)
        ch_versions = ch_versions.mix(CNVKIT_BATCH_TUMORONLY.out.versions.first())

    emit:
        versions = ch_versions         // channel: [ versions.yml ]

}
