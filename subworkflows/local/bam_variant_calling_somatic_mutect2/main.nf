//
// Run GATK mutect2 in tumor normal mode, getepileupsummaries, calculatecontamination, learnreadorientationmodel and filtermutectcalls
//

include { GATK4_CALCULATECONTAMINATION                                      } from '../../../modules/nf-core/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS                                           } from '../../../modules/nf-core/gatk4/filtermutectcalls/main'
include { GATK4_GATHERPILEUPSUMMARIES as GATK4_GATHERPILEUPSUMMARIES_NORMAL } from '../../../modules/nf-core/gatk4/gatherpileupsummaries/main'
include { GATK4_GATHERPILEUPSUMMARIES as GATK4_GATHERPILEUPSUMMARIES_TUMOR  } from '../../../modules/nf-core/gatk4/gatherpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES    as GATK4_GETPILEUPSUMMARIES_NORMAL    } from '../../../modules/nf-core/gatk4/getpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES    as GATK4_GETPILEUPSUMMARIES_TUMOR     } from '../../../modules/nf-core/gatk4/getpileupsummaries/main'
include { GATK4_LEARNREADORIENTATIONMODEL                                   } from '../../../modules/nf-core/gatk4/learnreadorientationmodel/main'
include { GATK4_MERGEMUTECTSTATS                                            } from '../../../modules/nf-core/gatk4/mergemutectstats/main'
include { GATK4_MERGEVCFS                                                   } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MUTECT2                                                     } from '../../../modules/nf-core/gatk4/mutect2/main'

workflow BAM_VARIANT_CALLING_SOMATIC_MUTECT2 {
    take:
    input                     // channel: [ meta, [ input ], [ input_index ] ]
    fasta                     // channel: /path/to/reference/fasta
    fai                       // channel: /path/to/reference/fasta/index
    dict                      // channel: /path/to/reference/fasta/dictionary
    germline_resource         // channel: /path/to/germline/resource
    germline_resource_tbi     // channel: /path/to/germline/index
    panel_of_normals          // channel: /path/to/panel/of/normals
    panel_of_normals_tbi      // channel: /path/to/panel/of/normals/index
    intervals                 // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals
    joint_mutect2             // boolean: [mandatory] [default: false] run mutect2 in joint mode

    main:
    versions = Channel.empty()

    germline_resource_pileup     = germline_resource_tbi ? germline_resource : Channel.empty()
    germline_resource_pileup_tbi = germline_resource_tbi ?: Channel.empty()

    // Combine input and intervals for spread and gather strategy
    input_intervals = input.combine(intervals)
        // Move num_intervals to meta map and reorganize channel for GATK4_MUTECT2 module
        .map{ meta, input_list, input_index_list, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], input_list, input_index_list, intervals ] }

    if (joint_mutect2) {
        // Separate normal cram files and remove duplicates
        ch_normal_cram = input.map{ meta, cram, crai -> [ meta - meta.subMap('tumor_id') + [id:meta.patient], cram[0], crai[0] ] }.unique()
        // Extract tumor cram files
        ch_tumor_cram = input.map{ meta, cram, crai -> [ meta - meta.subMap('tumor_id') + [id:meta.patient], cram[1], crai[1] ] }
        // Merge normal and tumor crams by patient
        ch_tn_cram = ch_normal_cram.mix(ch_tumor_cram).groupTuple()
        // Combine input and intervals for scatter and gather strategy
        ch_tn_intervals = ch_tn_cram.combine(intervals)
            // Move num_intervals to meta map and reorganize channel for GATK4_MUTECT2 module
            .map{ meta, cram, crai, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], cram, crai, intervals ] }
        GATK4_MUTECT2( ch_tn_intervals, fasta, fai, dict, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi)
    }
    else {
        // Perform variant calling using mutect2 module pair mode
        GATK4_MUTECT2( input_intervals, fasta, fai, dict, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi)
    }

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_branch = GATK4_MUTECT2.out.vcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more tbi(s) from the same sample
    tbi_branch = GATK4_MUTECT2.out.tbi.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more vcf(s) from the same sample
    stats_branch = GATK4_MUTECT2.out.stats.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more vcf(s) from the same sample
    f1r2_branch = GATK4_MUTECT2.out.f1r2.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    vcf_to_merge = vcf_branch.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ] }.groupTuple()
    stats_to_merge = stats_branch.intervals.map{ meta, stats -> [ groupKey(meta, meta.num_intervals), stats ] }.groupTuple()
    f1r2_to_merge = f1r2_branch.intervals.map{ meta, f1r2 -> [ groupKey(meta, meta.num_intervals), f1r2 ] }.groupTuple()

    GATK4_MERGEVCFS(vcf_to_merge, dict)
    GATK4_MERGEMUTECTSTATS(stats_to_merge)

    // Mix intervals and no_intervals channels together and remove no longer necessary field: normal_id, tumor_id, num_intervals
    vcf = Channel.empty().mix(GATK4_MERGEVCFS.out.vcf, vcf_branch.no_intervals).map{ meta, vcf -> [ meta - meta.subMap('normal_id', 'tumor_id', 'num_intervals'), vcf ]}
    tbi = Channel.empty().mix(GATK4_MERGEVCFS.out.tbi, tbi_branch.no_intervals).map{ meta, tbi -> [ meta - meta.subMap('normal_id', 'tumor_id', 'num_intervals'), tbi ]}
    stats = Channel.empty().mix(GATK4_MERGEMUTECTSTATS.out.stats, stats_branch.no_intervals).map{ meta, stats -> [ meta - meta.subMap('normal_id', 'tumor_id', 'num_intervals'), stats ]}
    f1r2 = Channel.empty().mix(f1r2_to_merge, f1r2_branch.no_intervals).map{ meta, f1r2 -> [ meta - meta.subMap('normal_id', 'tumor_id', 'num_intervals'), f1r2 ]}

    // Generate artifactpriors using learnreadorientationmodel on the f1r2 output of mutect2
    GATK4_LEARNREADORIENTATIONMODEL(f1r2)

    pileup = input_intervals.multiMap{  meta, input_list, input_index_list, intervals ->
        tumor: [ meta, input_list[1], input_index_list[1], intervals ]
        normal: [ meta, input_list[0], input_index_list[0], intervals ]
    }

    // Prepare input channel for normal pileup summaries.
    // Remember, the input channel contains tumor-normal pairs, so there will be multiple copies of the normal sample for each tumor for a given patient.
    // Therefore, we use unique function to generate normal pileup summaries once for each patient for better efficiency.
    pileup_normal = pileup.normal.map{ meta, cram, crai, intervals -> [ meta - meta.subMap('tumor_id') + [ id:meta.normal_id ], cram, crai, intervals] }.unique()
    // Prepare input channel for tumor pileup summaries.
    pileup_tumor = pileup.tumor.map{ meta, cram, crai, intervals -> [ meta + [ id:meta.tumor_id ], cram, crai, intervals ] }

    // Generate pileup summary tables using getepileupsummaries. tumor sample should always be passed in as the first input and input list entries of vcf_to_filter,
    GATK4_GETPILEUPSUMMARIES_NORMAL(pileup_normal, fasta, fai, dict, germline_resource_pileup, germline_resource_pileup_tbi)
    GATK4_GETPILEUPSUMMARIES_TUMOR(pileup_tumor, fasta, fai, dict, germline_resource_pileup, germline_resource_pileup_tbi)

    // Figuring out if there is one or more table(s) from the same sample
    pileup_table_normal_branch = GATK4_GETPILEUPSUMMARIES_NORMAL.out.table.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more table(s) from the same sample
    pileup_table_tumor_branch = GATK4_GETPILEUPSUMMARIES_TUMOR.out.table.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    pileup_table_normal_to_merge = pileup_table_normal_branch.intervals.map{ meta, table -> [ groupKey(meta, meta.num_intervals), table ] }.groupTuple()
    pileup_table_tumor_to_merge = pileup_table_tumor_branch.intervals.map{ meta, table -> [ groupKey(meta, meta.num_intervals), table ] }.groupTuple()

    // Merge Pileup Summaries
    GATK4_GATHERPILEUPSUMMARIES_NORMAL(pileup_table_normal_to_merge, dict.map{ meta, dict -> [ dict ] })
    GATK4_GATHERPILEUPSUMMARIES_TUMOR(pileup_table_tumor_to_merge, dict.map{ meta, dict -> [ dict ] })

    // Do some channel magic to generate tumor-normal pairs again.
    // This is necessary because we generated one normal pileup summary for each patient but we need run calculate contamination for each tumor-normal pair.
    pileup_table_tumor = Channel.empty().mix(GATK4_GATHERPILEUPSUMMARIES_TUMOR.out.table, pileup_table_tumor_branch.no_intervals).map{meta, table -> [ meta - meta.subMap('normal_id', 'tumor_id', 'num_intervals') + [id:meta.patient], meta.id, table ] }
    pileup_table_normal= Channel.empty().mix(GATK4_GATHERPILEUPSUMMARIES_NORMAL.out.table, pileup_table_normal_branch.no_intervals).map{meta, table -> [ meta - meta.subMap('normal_id', 'tumor_id', 'num_intervals') + [id:meta.patient], meta.id, table ] }
    ch_calculatecontamination_in_tables = pileup_table_tumor.combine(
        pileup_table_normal, by:0).map{
        meta, tumor_id, tumor_table, normal_id, normal_table -> [ meta + [ id: tumor_id + "_vs_" + normal_id ], tumor_table, normal_table]
        }

    GATK4_CALCULATECONTAMINATION(ch_calculatecontamination_in_tables)

    // Reduce the meta to only patient name if joint_mutect2 otherwise keep regular ID
    contamination_table = joint_mutect2 ? GATK4_CALCULATECONTAMINATION.out.contamination.map{ meta, cont -> [ meta - meta.subMap('tumor_id') + [id: meta.patient], cont]}.groupTuple() : GATK4_CALCULATECONTAMINATION.out.contamination
    segmentation_table  = joint_mutect2 ? GATK4_CALCULATECONTAMINATION.out.segmentation.map{  meta, seg  -> [ meta - meta.subMap('tumor_id') + [id: meta.patient], seg]}.groupTuple()  : GATK4_CALCULATECONTAMINATION.out.segmentation

    // Mutect2 calls filtered by filtermutectcalls using the artifactpriors, contamination and segmentation tables
    vcf_to_filter = vcf.join(tbi, failOnDuplicate: true, failOnMismatch: true)
        .join(stats, failOnDuplicate: true, failOnMismatch: true)
        .join(GATK4_LEARNREADORIENTATIONMODEL.out.artifactprior, failOnDuplicate: true, failOnMismatch: true)
        .join(segmentation_table, failOnDuplicate: true, failOnMismatch: true)
        .join(contamination_table, failOnDuplicate: true, failOnMismatch: true)
        .map{ meta, vcf, tbi, stats, orientation, seg, cont -> [ meta, vcf, tbi, stats, orientation, seg, cont, [] ] }

    GATK4_FILTERMUTECTCALLS(vcf_to_filter, fasta, fai, dict)

    vcf_filtered = GATK4_FILTERMUTECTCALLS.out.vcf
        // add variantcaller to meta map
        .map{ meta, vcf -> [ meta + [ variantcaller:'mutect2' ], vcf ] }

    versions = versions.mix(GATK4_CALCULATECONTAMINATION.out.versions)
    versions = versions.mix(GATK4_FILTERMUTECTCALLS.out.versions)
    versions = versions.mix(GATK4_GATHERPILEUPSUMMARIES_NORMAL.out.versions)
    versions = versions.mix(GATK4_GATHERPILEUPSUMMARIES_TUMOR.out.versions)
    versions = versions.mix(GATK4_GETPILEUPSUMMARIES_NORMAL.out.versions)
    versions = versions.mix(GATK4_GETPILEUPSUMMARIES_TUMOR.out.versions)
    versions = versions.mix(GATK4_LEARNREADORIENTATIONMODEL.out.versions)
    versions = versions.mix(GATK4_MERGEMUTECTSTATS.out.versions)
    versions = versions.mix(GATK4_MERGEVCFS.out.versions)
    versions = versions.mix(GATK4_MUTECT2.out.versions)

    emit:
    vcf   // channel: [ meta, vcf ]
    tbi   // channel: [ meta, tbi ]
    stats // channel: [ meta, stats ]

    vcf_filtered                                        // channel: [ meta, vcf ]
    index_filtered  = GATK4_FILTERMUTECTCALLS.out.tbi   // channel: [ meta, tbi ]
    stats_filtered  = GATK4_FILTERMUTECTCALLS.out.stats // channel: [ meta, stats ]

    artifact_priors = GATK4_LEARNREADORIENTATIONMODEL.out.artifactprior // channel: [ meta, artifactprior ]

    contamination_table // channel: [ meta, contamination ]
    pileup_table_normal // channel: [ meta, table_normal ]
    pileup_table_tumor  // channel: [ meta, table_tumor ]
    segmentation_table  // channel: [ meta, segmentation ]

    versions // channel: [ versions.yml ]
}
