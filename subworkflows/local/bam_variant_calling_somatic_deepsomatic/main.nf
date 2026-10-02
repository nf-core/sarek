//
// DEEPSOMATIC: tumour-normal somatic variant calling
//

include { DEEPSOMATIC } from '../../../modules/nf-core/deepsomatic/main'
include { GATK4_MERGEVCFS as MERGE_DEEPSOMATIC_VCF } from '../../../modules/nf-core/gatk4/mergevcfs/main'

workflow BAM_VARIANT_CALLING_SOMATIC_DEEPSOMATIC {
    take:
    cram // channel: [mandatory] [ meta, normal_cram, normal_crai, tumour_cram, tumour_crai ]
    dict // channel: [mandatory] [ meta, dict ]
    fasta // channel: [mandatory] [ meta, fasta ]
    fasta_fai // channel: [mandatory] [ meta, fasta_fai ]
    intervals // channel: [mandatory] [ interval, num_intervals ] or [ [], 0 ]

    main:
    cram_intervals = cram
        .combine(intervals)
        .map { meta, normal_cram, normal_crai, tumor_cram, tumor_crai, interval, num_intervals ->
            [meta + [num_intervals: num_intervals], normal_cram, normal_crai, tumor_cram, tumor_crai, interval]
        }

    deepsomatic_input = cram_intervals.multiMap { meta, normal_cram, normal_crai, tumor_cram, tumor_crai, interval ->
        cram: [meta, normal_cram, normal_crai, tumor_cram, tumor_crai]
        interval: [[id: interval ? interval.baseName : 'no_intervals'], interval]
    }

    DEEPSOMATIC(
        deepsomatic_input.cram,
        deepsomatic_input.interval,
        fasta,
        fasta_fai,
        [[id: 'no_gzi'], []],
    )

    vcf_branch = DEEPSOMATIC.out.vcf.branch { meta, _vcf ->
        intervals: meta.num_intervals > 1
        no_intervals: meta.num_intervals <= 1
    }

    tbi_branch = DEEPSOMATIC.out.vcf_tbi.branch { meta, _tbi ->
        intervals: meta.num_intervals > 1
        no_intervals: meta.num_intervals <= 1
    }

    vcf_to_merge = vcf_branch.intervals
        .map { meta, vcf -> [groupKey(meta, meta.num_intervals), vcf] }
        .groupTuple()

    MERGE_DEEPSOMATIC_VCF(vcf_to_merge, dict)

    vcf = channel.empty()
        .mix(MERGE_DEEPSOMATIC_VCF.out.vcf, vcf_branch.no_intervals)
        .map { meta, vcf -> [meta - meta.subMap('num_intervals') + [variantcaller: 'deepsomatic'], vcf] }

    tbi = channel.empty()
        .mix(MERGE_DEEPSOMATIC_VCF.out.tbi, tbi_branch.no_intervals)
        .map { meta, tbi -> [meta - meta.subMap('num_intervals') + [variantcaller: 'deepsomatic'], tbi] }

    emit:
    vcf // channel: [ meta, vcf.gz ]
    tbi // channel: [ meta, vcf.gz.tbi ]
}
