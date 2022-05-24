//
// CNV calling TUMOR_ONLY
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run


include {CNVKIT_ANTITARGET                          } from '../../../../../modules/nf-core/modules/cnvkit/antitarget/main'
include {CNVKIT_REFERENCE                           } from '../../../../../modules/nf-core/modules/cnvkit/reference/main'
include {CNVKIT_BATCH as CNVKIT_BATCH_TUMORONLY     } from '../../../../../modules/nf-core/modules/cnvkit/batch/main'

workflow RUN_CNVKIT_TUMORONLY {
    take:
        cram_recalibrated   // channel: [mandatory] cram tumor
        fasta               // channel: [mandatory] fasta
        targets             // channel: [mandatory] bed
        reference           // channel: [] cnn

    main:
        ch_versions = Channel.empty()

        // prepare a reference for tumor_only mode based on target_baits
        targets.view()

        CNVKIT_ANTITARGET(targets.map{ it -> [[id:it[0].baseName], it] })

        CNVKIT_REFERENCE(fasta, targets, CNVKIT_ANTITARGET.out.bed.map{ meta,bed -> [bed]} )

        // use reference for calling CNVs
        // cram_input needs the fasta reference genome for bam_conversion

        CNVKIT_BATCH_TUMORONLY(cram_recalibrated, fasta, [], CNVKIT_REFERENCE.out.cnn)

        ch_versions = ch_versions.mix(CNVKIT_ANTITARGET.out.versions)
        ch_versions = ch_versions.mix(CNVKIT_REFERENCE.out.versions)
        ch_versions = ch_versions.mix(CNVKIT_BATCH_TUMORONLY.out.versions)

    emit:
        versions = ch_versions         // channel: [ versions.yml ]

}
