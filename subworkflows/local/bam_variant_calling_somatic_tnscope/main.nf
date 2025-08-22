//
//
// SENTIEON TNSCOPE: tumor-normal mode variantcalling
//

include { SENTIEON_TNSCOPE                 } from '../../../modules/nf-core/sentieon/tnscope'
include { GATK4_MERGEVCFS as MERGE_TNSCOPE } from '../../../modules/nf-core/gatk4/mergevcfs'

workflow BAM_VARIANT_CALLING_SOMATIC_TNSCOPE {
    take:
    input                 // channel: [ meta, [ input ], [ input_index ] ]
    fasta                 // channel: [ meta, /path/to/reference/fasta ]
    fai                   // channel: [ meta, /path/to/reference/fasta/index ]
    dict                  // channel: [ meta, /path/to/reference/dict ]
    germline_resource     // channel: /path/to/germline/resource
    germline_resource_tbi // channel: /path/to/germline/index
    panel_of_normals      // channel: /path/to/panel/of/normals
    panel_of_normals_tbi  // channel: /path/to/panel/of/normals/index
    intervals             // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals

    main:
    versions = Channel.empty()

    // Combine input and intervals for spread and gather strategy
    // Move num_intervals to meta map and reorganize channel for TNSCOPE module
    input_intervals = input
        .combine(intervals)
        .map { meta, input_, index, intervals_, num_intervals -> [meta + [num_intervals: num_intervals], input_, index, intervals_] }

    SENTIEON_TNSCOPE(
        input_intervals,
        fasta,
        fai,
        germline_resource.map { resource -> [[id: "resource"], resource] },
        germline_resource_tbi.map { index -> [[id: "resource"], index] },
        panel_of_normals.map { pon -> [[id: "pon"], pon] },
        panel_of_normals_tbi.map { index -> [[id: "pon"], index] },
        [[], []],
        [[], []],
    )

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_branch = SENTIEON_TNSCOPE.out.vcf.branch {
        intervals: it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }
    // Figuring out if there is one or more tbi(s) from the same sample
    tbi_branch = SENTIEON_TNSCOPE.out.index.branch {
        intervals: it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    vcf_to_merge = vcf_branch.intervals.map { meta, vcf -> [groupKey(meta, meta.num_intervals), vcf] }.groupTuple()

    // Merge if required
    MERGE_TNSCOPE(vcf_to_merge, dict)

    // Mix intervals and no_intervals channels together
    // Remove unnecessary metadata and add variantcaller
    vcf = Channel.empty()
        .mix(MERGE_TNSCOPE.out.vcf, vcf_branch.no_intervals)
        .map { meta, vcf -> [meta - meta.subMap('num_intervals') + [variantcaller: 'sentieon_tnscope'], vcf] }

    index = Channel.empty()
        .mix(MERGE_TNSCOPE.out.tbi, tbi_branch.no_intervals)
        .map { meta, tbi -> [meta - meta.subMap('num_intervals') + [variantcaller: 'sentieon_tnscope'], tbi] }

    versions = versions.mix(SENTIEON_TNSCOPE.out.versions)
    versions = versions.mix(MERGE_TNSCOPE.out.versions)

    emit:
    vcf      // channel: [ meta, vcf ]
    index    // channel: [ meta, index ]
    versions // channel: [ versions.yml ]
}
