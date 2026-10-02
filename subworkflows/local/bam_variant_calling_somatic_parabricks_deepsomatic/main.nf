//
// PARABRICKS DEEPSOMATIC: GPU-accelerated tumour-normal somatic variant calling
//

include { PARABRICKS_DEEPSOMATIC } from '../../../modules/nf-core/parabricks/deepsomatic/main'

workflow BAM_VARIANT_CALLING_SOMATIC_PARABRICKS_DEEPSOMATIC {
    take:
    cram // channel: [mandatory] [ meta, normal_cram, normal_crai, tumour_cram, tumour_crai ]
    fasta // channel: [mandatory] [ meta, fasta ]
    intervals_bed_combined // channel: [optional] [] or [ intervals.bed ]

    main:
    intervals_bed = intervals_bed_combined.map { bed -> [bed ? bed[0] : []] }

    cram_intervals = cram
        .combine(intervals_bed)
        .map { meta, normal_cram, normal_crai, tumor_cram, tumor_crai, interval ->
            [meta, tumor_cram, tumor_crai, normal_cram, normal_crai, interval]
        }

    PARABRICKS_DEEPSOMATIC(cram_intervals, fasta)

    vcf = PARABRICKS_DEEPSOMATIC.out.vcf.map { meta, vcf -> [meta + [variantcaller: 'parabricks_deepsomatic'], vcf] }

    tbi = PARABRICKS_DEEPSOMATIC.out.vcf_tbi.map { meta, tbi -> [meta + [variantcaller: 'parabricks_deepsomatic'], tbi] }

    emit:
    vcf // channel: [ meta, vcf.gz ]
    tbi // channel: [ meta, vcf.gz.tbi ]
}
