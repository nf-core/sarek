
include { GATK4_CALCULATECONTAMINATION as CALCULATECONTAMINATION        } from '../../../modules/nf-core/gatk4/calculatecontamination/main'
include { GATK4_GATHERPILEUPSUMMARIES  as GATHERPILEUPSUMMARIES_NORMAL  } from '../../../modules/nf-core/gatk4/gatherpileupsummaries/main'
include { GATK4_GATHERPILEUPSUMMARIES  as GATHERPILEUPSUMMARIES_TUMOR   } from '../../../modules/nf-core/gatk4/gatherpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES     as GETPILEUPSUMMARIES_NORMAL     } from '../../../modules/nf-core/gatk4/getpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES     as GETPILEUPSUMMARIES_TUMOR      } from '../../../modules/nf-core/gatk4/getpileupsummaries/main'
// TO-DO: Remove the following out-commented include-statement - if it is not needed.
// include { GATK4_LEARNREADORIENTATIONMODEL as LEARNREADORIENTATIONMODEL    } from '../../../modules/nf-core/gatk4/learnreadorientationmodel/main'
include { GATK4_MERGEMUTECTSTATS       as MERGEMUTECTSTATS              } from '../../../modules/nf-core/gatk4/mergemutectstats/main'
include { GATK4_MERGEVCFS              as MERGE_TNHAPLOTYPER2           } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { SENTIEON_TNFILTER                                             } from '../../../modules/nf-core/sentieon/tnfilter/main'
include { SENTIEON_TNHAPLOTYPER2       as SENTIEON_TNHAPLOTYPER2_PAIRED } from '../../../modules/nf-core/sentieon/tnhaplotyper2/main'

workflow BAM_VARIANT_CALLING_SOMATIC_SENTIEON_TNHAPLOTYPER2 {
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
    joint_tnhaplotyper2       // boolean: [mandatory] [default: false] run mutect2 in joint mode

    main:
    versions = Channel.empty()

    //If no germline resource is provided, then create an empty channel to avoid GetPileupsummaries from being run
    germline_resource_pileup     = germline_resource_tbi ? germline_resource : Channel.empty()
    germline_resource_pileup_tbi = germline_resource_tbi ?: Channel.empty()

    // Combine input and intervals for spread and gather strategy
    input_intervals = input.combine(intervals)
        // Move num_intervals to meta map and reorganize channel for SENTIEON_TNHAPLOTYPER2 module
        .map{ meta, input_list, input_index_list, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], input_list, input_index_list, intervals ] }

    if (joint_tnhaplotyper2) {
        // Separate normal cram files
        // Extract tumor cram files
        ch_cram = input.multiMap{ meta, cram, crai ->
                normal: [ meta - meta.subMap('tumor_id') , cram[0], crai[0] ]
                tumor:  [ meta - meta.subMap('tumor_id') , cram[1], crai[1] ]
            }

        // Remove duplicates from normal channel and merge normal and tumor crams by patient
        ch_tn_cram =  ch_cram.normal.unique().mix(ch_cram.tumor).groupTuple()
        // Combine input and intervals for scatter and gather strategy
        ch_tn_intervals = ch_tn_cram.combine(intervals)
            // Move num_intervals to meta map and reorganize channel for SENTIEON_TNHAPLOTYPER2 module
            // meta: [id:patient_id, num_intervals, patient, sex]
            .map{ meta, cram, crai, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], cram, crai, intervals ] }

        SENTIEON_TNHAPLOTYPER2_PAIRED(
            ch_tn_intervals,
            dict,
            fasta,
            fai,
            germline_resource.map{ it -> [ [ id:it.baseName ], it ] },
            germline_resource_tbi.map{ it -> [ [ id:it.baseName ], it ] },
            panel_of_normals.map{ it -> [ [ id:it.baseName ], it ] },
            panel_of_normals_tbi.map{ it -> [ [ id:it.baseName ], it ] },
            true,  // TO-DO: These things shouldn't be hardcoded
            true
        )
    } else {
        // Perform variant calling using mutect2 module pair mode
        // meta: [id:tumor_id_vs_normal_id, normal_id, num_intervals, patient, sex, tumor_id]
        SENTIEON_TNHAPLOTYPER2_PAIRED(
            input_intervals,
            dict,
            fasta,
            fai,
            germline_resource.map{ it -> [ [ id:it.baseName ], it ] },
            germline_resource_tbi.map{ it -> [ [ id:it.baseName ], it ] },
            panel_of_normals.map{ it -> [ [ id:it.baseName ], it ] },
            panel_of_normals_tbi.map{ it -> [ [ id:it.baseName ], it ] },
            true,  // TO-DO: These things shouldn't be hardcoded
            true
        )
    }

    vcf_to_filter = SENTIEON_TNHAPLOTYPER2_PAIRED.out.vcf
        .join(SENTIEON_TNHAPLOTYPER2_PAIRED.out.index, failOnDuplicate: true, failOnMismatch: true)
        .join(SENTIEON_TNHAPLOTYPER2_PAIRED.out.stats, failOnDuplicate: true, failOnMismatch: true)
        .join(SENTIEON_TNHAPLOTYPER2_PAIRED.out.contamination_data, failOnDuplicate: true, failOnMismatch: true)
        .join(SENTIEON_TNHAPLOTYPER2_PAIRED.out.contamination_segments, failOnDuplicate: true, failOnMismatch: true)
        .join(SENTIEON_TNHAPLOTYPER2_PAIRED.out.orientation_data, failOnDuplicate: true, failOnMismatch: true)
        // .map{ meta, vcf, tbi, stats, seg, cont -> [ meta, vcf, tbi, stats, [], seg, cont, [] ] }

/*
    TO-DO: Clean up the following
    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_branch = SENTIEON_TNHAPLOTYPER2_PAIRED.out.vcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more tbi(s) from the same sample
    tbi_branch = SENTIEON_TNHAPLOTYPER2_PAIRED.out.index.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more vcf(s) from the same sample
    stats_branch = SENTIEON_TNHAPLOTYPER2_PAIRED.out.stats.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    /* TO-DO: Doesn't seem relevant for tnhaplotyper2
    // Figuring out if there is one or more vcf(s) from the same sample
    f1r2_branch = SENTIEON_TNHAPLOTYPER2.out.f1r2.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }
    */

/*
    // Only when using intervals
    vcf_to_merge = vcf_branch.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ] }.groupTuple()
    stats_to_merge = stats_branch.intervals.map{ meta, stats -> [ groupKey(meta, meta.num_intervals), stats ] }.groupTuple()
    // TO-DO: Doesn't seem relevant for tnhaplotyper2
    // f1r2_to_merge = f1r2_branch.intervals.map{ meta, f1r2 -> [ groupKey(meta, meta.num_intervals), f1r2 ] }.groupTuple()

    MERGE_TNHAPLOTYPER2(vcf_to_merge, dict)
    MERGEMUTECTSTATS(stats_to_merge)

    // Mix intervals and no_intervals channels together and remove no longer necessary field: normal_id, tumor_id, num_intervals
    vcf = Channel.empty().mix(MERGE_TNHAPLOTYPER2.out.vcf, vcf_branch.no_intervals).map{ meta, vcf ->
            [ joint_mutect2 ? meta - meta.subMap('normal_id', 'num_intervals') :  meta - meta.subMap('num_intervals') , vcf ]
    }
    tbi = Channel.empty().mix(MERGE_TNHAPLOTYPER2.out.tbi, tbi_branch.no_intervals).map{ meta, tbi->
            [ joint_mutect2 ? meta - meta.subMap('normal_id', 'num_intervals') :  meta - meta.subMap('num_intervals'), tbi ]
    }
    stats = Channel.empty().mix(MERGEMUTECTSTATS.out.stats, stats_branch.no_intervals).map{ meta, stats ->
            [ joint_mutect2 ? meta - meta.subMap('normal_id', 'num_intervals') :  meta - meta.subMap('num_intervals'), stats ]
    }
    // TO-DO: Doesn't seem relevant for tnhaplotyper2
    /*
    f1r2 = Channel.empty().mix(f1r2_to_merge, f1r2_branch.no_intervals).map{ meta, f1r2->
            [ joint_mutect2 ? meta - meta.subMap('normal_id', 'num_intervals') :  meta - meta.subMap('num_intervals') , f1r2 ]
    }
    */

    // TO-DO: Doesn't seem relevant for tnhaplotyper2
    // Generate artifactpriors using learnreadorientationmodel on the f1r2 output of mutect2
    // LEARNREADORIENTATIONMODEL(f1r2)
/*
    pileup = input_intervals.multiMap{  meta, input_list, input_index_list, intervals ->
        tumor: [ meta, input_list[1], input_index_list[1], intervals ]
        normal: [ meta, input_list[0], input_index_list[0], intervals ]
    }

    // Prepare input channel for normal pileup summaries.
    // Remember, the input channel contains tumor-normal pairs, so there will be multiple copies of the normal sample for each tumor for a given patient.
    // Therefore, we use unique function to generate normal pileup summaries once for each patient for better efficiency.
    pileup_normal = pileup.normal.map{ meta, cram, crai, intervals -> [ meta - meta.subMap('tumor_id') + [ id:meta.normal_id ], cram, crai, intervals] }.unique()
    // Prepare input channel for tumor pileup summaries.
    pileup_tumor = pileup.tumor.map{ meta, cram, crai, intervals -> [ meta - meta.subMap('normal_id') + [ id:meta.tumor_id ], cram, crai, intervals ] }

    // Generate pileup summary tables using getepileupsummaries. tumor sample should always be passed in as the first input and input list entries of vcf_to_filter,
    GETPILEUPSUMMARIES_NORMAL(pileup_normal, fasta, fai, dict, germline_resource_pileup, germline_resource_pileup_tbi)
    GETPILEUPSUMMARIES_TUMOR(pileup_tumor, fasta, fai, dict, germline_resource_pileup, germline_resource_pileup_tbi)

    // Figuring out if there is one or more table(s) from the same sample
    pileup_table_normal_branch = GETPILEUPSUMMARIES_NORMAL.out.table.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more table(s) from the same sample
    pileup_table_tumor_branch = GETPILEUPSUMMARIES_TUMOR.out.table.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    pileup_table_normal_to_merge = pileup_table_normal_branch.intervals.map{ meta, table -> [ groupKey(meta, meta.num_intervals), table ] }.groupTuple()
    pileup_table_tumor_to_merge = pileup_table_tumor_branch.intervals.map{ meta, table -> [ groupKey(meta, meta.num_intervals), table ] }.groupTuple()

    // Merge Pileup Summaries
    GATHERPILEUPSUMMARIES_NORMAL(pileup_table_normal_to_merge, dict.map{ meta, dict -> [ dict ] })
    GATHERPILEUPSUMMARIES_TUMOR(pileup_table_tumor_to_merge, dict.map{ meta, dict -> [ dict ] })

    // Do some channel magic to generate tumor-normal pairs again.
    // This is necessary because we generated one normal pileup summary for each patient but we need run calculate contamination for each tumor-normal pair.
    pileup_table_tumor  = Channel.empty().mix(GATHERPILEUPSUMMARIES_TUMOR.out.table, pileup_table_tumor_branch.no_intervals).map{meta, table -> [ meta - meta.subMap('normal_id', 'tumor_id', 'num_intervals') + [id:meta.patient], meta.id, table ] }
    pileup_table_normal = Channel.empty().mix(GATHERPILEUPSUMMARIES_NORMAL.out.table, pileup_table_normal_branch.no_intervals).map{meta, table -> [ meta - meta.subMap('normal_id', 'tumor_id', 'num_intervals') + [id:meta.patient], meta.id, table ] }

    ch_calculatecontamination_in_tables = pileup_table_tumor.combine(
        pileup_table_normal, by:0).map{
        meta, tumor_id, tumor_table, normal_id, normal_table ->
            if(joint_mutect2){
                [ meta + [ id: tumor_id + "_vs_" + normal_id], tumor_table, normal_table]
            }else{
                // we need tumor and normal ID for further post processing
                [ meta + [ id: tumor_id + "_vs_" + normal_id, normal_id:normal_id, tumor_id:tumor_id ], tumor_table, normal_table]
            }
        }

    CALCULATECONTAMINATION(ch_calculatecontamination_in_tables)
*/
    // Initialize empty channel: Contamination calculation is run on pileup table, pileup is not run if germline resource is not provided
    calculatecontamination_out_seg = Channel.empty()
    calculatecontamination_out_cont = Channel.empty()
/*
    if (joint_mutect2) {
        // Reduce the meta to only patient name
        calculatecontamination_out_seg = CALCULATECONTAMINATION.out.segmentation.map{ meta, seg -> [ meta + [id: meta.patient], seg]}.groupTuple()
        calculatecontamination_out_cont = CALCULATECONTAMINATION.out.contamination.map{ meta, cont -> [ meta + [id: meta.patient], cont]}.groupTuple()
    }
    else {
        // Keep tumor_vs_normal ID
        calculatecontamination_out_seg = CALCULATECONTAMINATION.out.segmentation
        calculatecontamination_out_cont = CALCULATECONTAMINATION.out.contamination
    }

    // Mutect2 calls filtered by filtermutectcalls using the artifactpriors, contamination and segmentation tables
    // meta joint calling:  [id:patient_id, patient, sex]
    // meta paired calling: [id:tumorID_vs_normalID, normal_ID, patient, sex, tumorID]
    vcf_to_filter = vcf.join(tbi, failOnDuplicate: true, failOnMismatch: true)
        .join(stats, failOnDuplicate: true, failOnMismatch: true)
        // .join(LEARNREADORIENTATIONMODEL.out.artifactprior, failOnDuplicate: true, failOnMismatch: true)
        .join(calculatecontamination_out_seg)
        .join(calculatecontamination_out_cont)
        // .map{ meta, vcf, tbi, stats, orientation, seg, cont -> [ meta, vcf, tbi, stats, orientation, seg, cont, [] ] }
        .map{ meta, vcf, tbi, stats, seg, cont -> [ meta, vcf, tbi, stats, [], seg, cont, [] ] }

*/

    SENTIEON_TNFILTER(vcf_to_filter, fasta, fai)

    vcf_filtered = SENTIEON_TNFILTER.out.vcf
        // add variantcaller to meta map
        .map{ meta, vcf -> [ meta + [ variantcaller:'sentieon_tnhaplotyper2' ], vcf ] }

    /* INFO: Some of the following may be needed later...
       TO-DO: Clean up the following
    versions = versions.mix(MERGE_TNHAPLOTYPER2.out.versions)
    versions = versions.mix(CALCULATECONTAMINATION.out.versions)
    */
    versions = versions.mix(SENTIEON_TNFILTER.out.versions)
    /*
    versions = versions.mix(GETPILEUPSUMMARIES_NORMAL.out.versions)
    versions = versions.mix(GETPILEUPSUMMARIES_TUMOR.out.versions)
    versions = versions.mix(GATHERPILEUPSUMMARIES_NORMAL.out.versions)
    versions = versions.mix(GATHERPILEUPSUMMARIES_TUMOR.out.versions)
    */
    // TO-DO: Doesn't seem relevant for tnhaplotyper2
    // versions = versions.mix(LEARNREADORIENTATIONMODEL.out.versions)
    // INFO: May be needed later...
    // versions = versions.mix(MERGEMUTECTSTATS.out.versions)
    versions = versions.mix(SENTIEON_TNHAPLOTYPER2_PAIRED.out.versions)

    emit:
    // INFO: vcf and stats before came from vcf = Channel.empty().mix(MERGE_TNHAPLOTYPER2.out.vcf etc.
    vcf = SENTIEON_TNHAPLOTYPER2_PAIRED.out.vcf  // channel: [ meta, vcf ]
    stats = SENTIEON_TNHAPLOTYPER2_PAIRED.out.stats // channel: [ meta, stats ]

    vcf_filtered                                    // channel: [ meta, vcf ]
    index_filtered = SENTIEON_TNFILTER.out.vcf_tbi  // channel: [ meta, tbi ]
    stats_filtered = SENTIEON_TNFILTER.out.stats    // channel: [ meta, stats ]

    // TO-DO: Doesn't seem relevant for tnhaplotyper2
    // artifact_priors        = LEARNREADORIENTATIONMODEL.out.artifactprior // channel: [ meta, artifactprior ]

    // INFO: Temporarily just set to empty channels
    pileup_table_normal = Channel.empty() // channel: [ meta, table_normal ]
    pileup_table_tumor  = Channel.empty() // channel: [ meta, table_tumor ]

    contamination_table    = calculatecontamination_out_cont    // channel: [ meta, contamination ]
    segmentation_table     = calculatecontamination_out_seg     // channel: [ meta, segmentation ]

    versions // channel: [ versions.yml ]
}
