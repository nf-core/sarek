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
    bam_tumor // channel: [mandatory] [ meta, tumor_bam, tumor_bai]
    fasta // channel: [mandatory] [ meta, fasta ]
    dbsnp // channel: [mandatory] [ dbsnp ]
    dbsnp_tbi // channel: [mandatory] [ dbsnp_tbi ]

    main:
    // Combine normal and tumor data
    ch_bam = bam_tumor.join(bam_normal, by: [0])

    ch_call_input = ch_bam
        .combine(fasta)
        .map { meta_tumor_bam, tumor_bam, tumor_bai, _meta_normal_bam, normal_bam, normal_bai, _meta_fasta, fasta_file ->
            [meta_tumor_bam, tumor_bam, tumor_bai, normal_bam, normal_bai, fasta_file]
        }

    MUSE_CALL(ch_call_input)

    ch_sump_input = MUSE_CALL.out.txt.combine(dbsnp.combine(dbsnp_tbi))

    MUSE_SUMP(ch_sump_input)

    // add variantcaller to meta map
    vcf = MUSE_SUMP.out.vcf.map { meta, vcf -> [meta + [variantcaller: 'muse'], vcf] }
    tbi = MUSE_SUMP.out.tbi.map { meta, tbi -> [meta + [variantcaller: 'muse'], tbi] }

    emit:
    vcf
    tbi
}
