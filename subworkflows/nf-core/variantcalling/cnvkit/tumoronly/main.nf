//
// CNV calling TUMOR_ONLY
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include {CNVKIT_BATCH as CNVKIT_BATCH_TUMORONLY     } from '../../../../../modules/nf-core/modules/cnvkit/batch/main'

workflow RUN_CNVKIT_TUMORONLY {
    take:
        cram_recalibrated   // channel: [mandatory] cram tumor
        fasta               // channel: [mandatory] fasta
        fasta_fai           // channel: [optional]  fasta_fai
        targets             // channel: [mandatory] bed
        reference           // channel: [] cnn

    main:
        ch_versions = Channel.empty()

        // use reference for calling CNVs
        // cram_input needs the fasta reference genome for bam_conversion

        CNVKIT_BATCH_TUMORONLY(cram_recalibrated, fasta, fasta_fai, [], reference)

        ch_versions = ch_versions.mix(CNVKIT_BATCH_TUMORONLY.out.versions)

    emit:
        versions = ch_versions         // channel: [ versions.yml ]

}
