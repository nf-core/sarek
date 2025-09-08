//
// MuSE tumor-normal variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { MUSE_CALL } from '../../../modules/nf-core/muse/call'
include { MUSE_SUMP } from '../../../modules/nf-core/muse/sump'

workflow BAM_VARIANT_CALLING_SOMATIC_MUSE {
    take:
    bam_normal // channel: [mandatory] [ meta, normal_bam, normal_bai]
    bam_tumor  // channel: [mandatory] [ meta, tumor_bam, tumor_bai]
    fasta      // channel: [mandatory] [ meta, fasta ]
    dbsnp      // channel: [mandatory] [ dbsnp ]
    dbsnp_tbi  // channel: [mandatory] [ dbsnp_tbi ]

    main:
    versions = Channel.empty()
    ch_dbsnp = dbsnp.map { it -> [[id: it.name], it] }
    ch_dbsnp_tbi = dbsnp_tbi.map { it -> [[id: it.baseName], it] }
    ch_dbsnp_with_tbi = ch_dbsnp.join(ch_dbsnp_tbi, by: [0]).collect()
    // Join by meta.id

    // Combine normal and tumor data
    ch_bam = bam_tumor.join(bam_normal, by: [0])
    // Join by meta

    MUSE_CALL(
        ch_bam,
        fasta,
    )

    MUSE_SUMP(
        MUSE_CALL.out.txt,
        ch_dbsnp_with_tbi,
    )

    // add variantcaller to meta map
    vcf = MUSE_SUMP.out.vcf.map { meta, vcf -> [meta + [variantcaller: 'muse'], vcf] }

    versions = versions.mix(MUSE_CALL.out.versions)
    versions = versions.mix(MUSE_SUMP.out.versions)

    emit:
    vcf
    versions
}
