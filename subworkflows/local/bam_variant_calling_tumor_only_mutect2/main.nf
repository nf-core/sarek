//
// GATK MUTECT2 in tumor only mode, filtering in BAM_VCF_FILTER_MUTECT_CALLS_GATK
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BAM_VCF_FILTER_MUTECT_CALLS_GATK           } from '../bam_vcf_filter_mutect_calls_gatk/main'
include { GATK4_MERGEMUTECTSTATS as MERGEMUTECTSTATS } from '../../../modules/nf-core/gatk4/mergemutectstats'
include { GATK4_MERGEVCFS as MERGE_MUTECT2           } from '../../../modules/nf-core/gatk4/mergevcfs'
include { GATK4_MUTECT2 as MUTECT2                   } from '../../../modules/nf-core/gatk4/mutect2'

workflow BAM_VARIANT_CALLING_TUMOR_ONLY_MUTECT2 {
    take:
    input                      // channel: [ meta, [ input ], [ input_index ] ]
    fasta                      // channel: /path/to/reference/fasta
    fai                        // channel: /path/to/reference/fasta/index
    dict                       // channel: /path/to/reference/fasta/dictionary
    germline_resource          // channel: /path/to/germline/resource
    germline_resource_tbi      // channel: /path/to/germline/index
    panel_of_normals           // channel: /path/to/panel/of/normals
    panel_of_normals_tbi       // channel: /path/to/panel/of/normals/index
    intervals                  // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals
    joint_mutect2              // boolean: [mandatory] [default: false] run mutect2 in joint mode
    germline_resource_provided // boolean: [mandatory] germline resource provided, run contamination estimation

    main:

    // Combine input and intervals for spread and gather strategy
    input_intervals = input
        .combine(intervals)
        .map { meta, input_, index, intervals_, num_intervals -> [meta + [num_intervals: num_intervals], input_, index, intervals_] }

    if (joint_mutect2) {
        // Perform variant calling using mutect2 module in tumor single mode
        // Group cram files by patient
        input_joint = input
            .map { meta, input_, index -> [meta - meta.subMap('sample') + [id: meta.patient], input_, index] }
            .groupTuple()

        // Add intervals for scatter-gather scaling
        input_joint_intervals = input_joint
            .combine(intervals)
            .map { meta, cram, crai, intervals_, num_intervals -> [meta + [num_intervals: num_intervals], cram, crai, intervals_] }
        MUTECT2(input_joint_intervals, fasta, fai.map { meta, index -> [ meta, index, [] ] }, dict, [], [], germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi)
    }
    else {
        // Perform variant calling using mutect2 module in tumor single mode
        MUTECT2(input_intervals, fasta, fai.map { meta, index -> [ meta, index, [] ] }, dict, [], [], germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi)
    }

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_branch = MUTECT2.out.vcf.branch { meta, _vcf ->
        intervals: meta.num_intervals > 1
        no_intervals: meta.num_intervals <= 1
    }

    // Figuring out if there is one or more tbi(s) from the same sample
    tbi_branch = MUTECT2.out.tbi.branch { meta, _tbi ->
        intervals: meta.num_intervals > 1
        no_intervals: meta.num_intervals <= 1
    }

    // Figuring out if there is one or more stats(s) from the same sample
    stats_branch = MUTECT2.out.stats.branch { meta, _stats ->
        intervals: meta.num_intervals > 1
        no_intervals: meta.num_intervals <= 1
    }

    // Figuring out if there is one or more f1r2(s) from the same sample
    f1r2_branch = MUTECT2.out.f1r2.branch { meta, _f1r2 ->
        intervals: meta.num_intervals > 1
        no_intervals: meta.num_intervals <= 1
    }

    // Only when using intervals
    vcf_to_merge = vcf_branch.intervals.map { meta, vcf -> [groupKey(meta, meta.num_intervals), vcf] }.groupTuple()
    stats_to_merge = stats_branch.intervals.map { meta, stats -> [groupKey(meta, meta.num_intervals), stats] }.groupTuple()
    f1r2_to_merge = f1r2_branch.intervals.map { meta, f1r2 -> [groupKey(meta, meta.num_intervals), f1r2] }.groupTuple()

    MERGE_MUTECT2(vcf_to_merge, dict)
    MERGEMUTECTSTATS(stats_to_merge)

    // Mix intervals and no_intervals channels together
    // Remove unnecessary metadata
    vcf = channel.empty().mix(MERGE_MUTECT2.out.vcf, vcf_branch.no_intervals).map { meta, vcf -> [meta - meta.subMap('num_intervals'), vcf] }
    tbi = channel.empty().mix(MERGE_MUTECT2.out.tbi, tbi_branch.no_intervals).map { meta, tbi -> [meta - meta.subMap('num_intervals'), tbi] }
    stats = channel.empty().mix(MERGEMUTECTSTATS.out.stats, stats_branch.no_intervals).map { meta, stats -> [meta - meta.subMap('num_intervals'), stats] }
    f1r2 = channel.empty().mix(f1r2_to_merge, f1r2_branch.no_intervals).map { meta, f1r2 -> [meta - meta.subMap('num_intervals'), f1r2] }

    // Mutect2 calls filtered by filtermutectcalls in the shared filtering subworkflow
    calls = vcf
        .join(tbi, failOnDuplicate: true, failOnMismatch: true)
        .join(stats, failOnDuplicate: true, failOnMismatch: true)
        .join(f1r2, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, vcf_, tbi_, stats_, f1r2_ -> [meta + [variantcaller: 'mutect2'], vcf_, tbi_, stats_, f1r2_] }

    // One row per tumor sample, also in joint mode
    cram = input.map { meta, cram_, crai_ -> [meta + [variantcaller: 'mutect2'], [], [], cram_, crai_] }

    BAM_VCF_FILTER_MUTECT_CALLS_GATK(
        calls,
        cram,
        intervals,
        fasta,
        fai,
        dict,
        germline_resource,
        germline_resource_tbi,
        germline_resource_provided,
        joint_mutect2,
    )

    emit:
    vcf                 = BAM_VCF_FILTER_MUTECT_CALLS_GATK.out.vcf                 // channel: [ meta, vcf ]
    tbi                 = BAM_VCF_FILTER_MUTECT_CALLS_GATK.out.tbi                 // channel: [ meta, tbi ]
    stats_filtered      = BAM_VCF_FILTER_MUTECT_CALLS_GATK.out.stats_filtered      // channel: [ meta, stats ]
    artifact_priors     = BAM_VCF_FILTER_MUTECT_CALLS_GATK.out.artifact_priors     // channel: [ meta, artifactprior ]
    pileup_table        = BAM_VCF_FILTER_MUTECT_CALLS_GATK.out.pileup_table_tumor  // channel: [ meta, table ]
    contamination_table = BAM_VCF_FILTER_MUTECT_CALLS_GATK.out.contamination_table // channel: [ meta, contamination ]
    segmentation_table  = BAM_VCF_FILTER_MUTECT_CALLS_GATK.out.segmentation_table  // channel: [ meta, segmentation ]
}
