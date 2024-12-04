//
// MuSE tumor-normal variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { MUSE_CALL } from '../../../modules/nf-core/muse/call/main'
include { MUSE_SUMP } from '../../../modules/nf-core/muse/sump/main'

workflow BAM_VARIANT_CALLING_SOMATIC_MUSE {
    take:
    cram          // channel: [mandatory] [ meta, normal_cram, normal_crai, tumor_cram, tumor_crai ]
    fasta         // channel: [mandatory] [ meta, fasta ]
    dbsnp         // channel: [optional] [ dbsnp ]
    dbsnp_tbi     // channel: [optional] [ dbsnp_tbi ]

    main:
    versions = Channel.empty()

    MUSE_CALL(
        cram,
        fasta
    )

    MUSE_SUMP(
        MUSE_CALL.out.txt,
        dbsnp.map{ it -> [ [ id:it.baseName ], it, dbsnp_tbi ] }
    )

    // Mix intervals and no_intervals channels together
    vcf = Channel.empty().mix(MUSE_SUMP.out.vcf)
        .map{ meta, vcf -> [ meta + [ variantcaller: 'muse' ], vcf ] }

    versions = versions.mix(MUSE_CALL.out.versions)
    versions = versions.mix(MUSE_SUMP.out.versions)

    emit:
    vcf
    versions
}
