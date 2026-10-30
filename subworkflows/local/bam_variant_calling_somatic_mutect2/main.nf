//
// MUTECT2: tumour-normal mode variant calling, filtering in BAM_VCF_FILTER_MUTECT_CALLS_GATK
//

include { BAM_VCF_FILTER_MUTECT_CALLS_GATK           } from '../bam_vcf_filter_mutect_calls_gatk/main'
include { GATK4_MERGEMUTECTSTATS as MERGEMUTECTSTATS } from '../../../modules/nf-core/gatk4/mergemutectstats'
include { GATK4_MERGEVCFS as MERGE_MUTECT2           } from '../../../modules/nf-core/gatk4/mergevcfs'
include { GATK4_MUTECT2 as MUTECT2_PAIRED            } from '../../../modules/nf-core/gatk4/mutect2'

workflow BAM_VARIANT_CALLING_SOMATIC_MUTECT2 {
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
    //   Move num_intervals to meta map and reorganize channel for MUTECT2_PAIRED module
    input_intervals = input
        .combine(intervals)
        .map { meta, input_list, input_index_list, intervals_, num_intervals -> [meta + [num_intervals: num_intervals], input_list, input_index_list, intervals_] }

    if (joint_mutect2) {

        // Separate normal cram files
        // Extract tumor cram files
        ch_cram = input.multiMap { meta, cram, crai ->
            normal: [meta - meta.subMap('tumor_id'), cram[0], crai[0]]
            tumor: [meta - meta.subMap('tumor_id'), cram[1], crai[1]]
        }

        // Remove duplicates from normal channel and merge normal and tumor crams by patient
        ch_tn_cram = ch_cram.normal.unique().mix(ch_cram.tumor).groupTuple()
        // Combine input and intervals for scatter and gather strategy
        //   Move num_intervals to meta map and reorganize channel for MUTECT2_PAIRED module
        ch_tn_intervals = ch_tn_cram
            .combine(intervals)
            .map { meta, cram, crai, intervals_, num_intervals -> [meta + [num_intervals: num_intervals], cram, crai, intervals_] }

        MUTECT2_PAIRED(ch_tn_intervals, fasta, fai.map { meta, index -> [ meta, index, [] ] }, dict, [], [], germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi)
    }
    else {

        // Perform variant calling using mutect2 module pair mode
        // meta: [id:tumor_id_vs_normal_id, normal_id, num_intervals, patient, sex, tumor_id]
        MUTECT2_PAIRED(input_intervals, fasta, fai.map { meta, index -> [ meta, index, [] ] }, dict, [], [], germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi)
    }

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_branch = MUTECT2_PAIRED.out.vcf.branch { meta, _vcf ->
        intervals: meta.num_intervals > 1
        no_intervals: meta.num_intervals <= 1
    }

    // Figuring out if there is one or more tbi(s) from the same sample
    tbi_branch = MUTECT2_PAIRED.out.tbi.branch { meta, _tbi ->
        intervals: meta.num_intervals > 1
        no_intervals: meta.num_intervals <= 1
    }

    // Figuring out if there is one or more vcf(s) from the same sample
    stats_branch = MUTECT2_PAIRED.out.stats.branch { meta, _stats ->
        intervals: meta.num_intervals > 1
        no_intervals: meta.num_intervals <= 1
    }

    // Figuring out if there is one or more vcf(s) from the same sample
    f1r2_branch = MUTECT2_PAIRED.out.f1r2.branch { meta, _f1r2 ->
        intervals: meta.num_intervals > 1
        no_intervals: meta.num_intervals <= 1
    }

    // Only when using intervals
    vcf_to_merge = vcf_branch.intervals.map { meta, vcf -> [groupKey(meta, meta.num_intervals), vcf] }.groupTuple()
    stats_to_merge = stats_branch.intervals.map { meta, stats -> [groupKey(meta, meta.num_intervals), stats] }.groupTuple()
    f1r2_to_merge = f1r2_branch.intervals.map { meta, f1r2 -> [groupKey(meta, meta.num_intervals), f1r2] }.groupTuple()

    MERGE_MUTECT2(vcf_to_merge, dict)
    MERGEMUTECTSTATS(stats_to_merge)

    // Mix intervals and no_intervals channels together and remove no longer necessary field: normal_id, tumor_id, num_intervals
    vcf = channel.empty()
        .mix(MERGE_MUTECT2.out.vcf, vcf_branch.no_intervals)
        .map { meta, vcf ->
            [joint_mutect2 ? meta - meta.subMap('normal_id', 'num_intervals') : meta - meta.subMap('num_intervals'), vcf]
        }
    tbi = channel.empty()
        .mix(MERGE_MUTECT2.out.tbi, tbi_branch.no_intervals)
        .map { meta, tbi ->
            [joint_mutect2 ? meta - meta.subMap('normal_id', 'num_intervals') : meta - meta.subMap('num_intervals'), tbi]
        }
    stats = channel.empty()
        .mix(MERGEMUTECTSTATS.out.stats, stats_branch.no_intervals)
        .map { meta, stats ->
            [joint_mutect2 ? meta - meta.subMap('normal_id', 'num_intervals') : meta - meta.subMap('num_intervals'), stats]
        }
    f1r2 = channel.empty()
        .mix(f1r2_to_merge, f1r2_branch.no_intervals)
        .map { meta, f1r2 ->
            [joint_mutect2 ? meta - meta.subMap('normal_id', 'num_intervals') : meta - meta.subMap('num_intervals'), f1r2]
        }

    // Mutect2 calls filtered by filtermutectcalls in the shared filtering subworkflow
    // meta joint calling:  [id:patient_id, patient, sex, variantcaller]
    // meta paired calling: [id:tumorID_vs_normalID, normal_ID, patient, sex, tumorID, variantcaller]
    calls = vcf
        .join(tbi, failOnDuplicate: true, failOnMismatch: true)
        .join(stats, failOnDuplicate: true, failOnMismatch: true)
        .join(f1r2, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, vcf_, tbi_, stats_, f1r2_ -> [meta + [variantcaller: 'mutect2'], vcf_, tbi_, stats_, f1r2_] }

    // One row per tumour/normal pair, also in joint mode
    cram = input.map { meta, input_list, input_index_list -> [meta + [variantcaller: 'mutect2'], input_list[0], input_index_list[0], input_list[1], input_index_list[1]] }

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
    pileup_table_normal = BAM_VCF_FILTER_MUTECT_CALLS_GATK.out.pileup_table_normal // channel: [ meta, table_normal ]
    pileup_table_tumor  = BAM_VCF_FILTER_MUTECT_CALLS_GATK.out.pileup_table_tumor  // channel: [ meta, table_tumor ]
    contamination_table = BAM_VCF_FILTER_MUTECT_CALLS_GATK.out.contamination_table // channel: [ meta, contamination ]
    segmentation_table  = BAM_VCF_FILTER_MUTECT_CALLS_GATK.out.segmentation_table  // channel: [ meta, segmentation ]
}
