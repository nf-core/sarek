//
// Run GATK mutect2 in multi sample mode, getpileupsummaries, calculatecontamination, learnreadorientationmodel and filtermutectcalls
//

include { GATK4_MUTECT2                   as MUTECT2_PAIRED                } from '../../../modules/nf-core/gatk4/mutect2/main'
include { GATK4_MERGEVCFS                 as MERGE_MUTECT2                 } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEMUTECTSTATS          as MERGEMUTECTSTATS              } from '../../../modules/nf-core/gatk4/mergemutectstats/main'
include { GATK4_LEARNREADORIENTATIONMODEL as LEARNREADORIENTATIONMODEL     } from '../../../modules/nf-core/gatk4/learnreadorientationmodel/main'
include { GATK4_GATHERPILEUPSUMMARIES     as GATHERPILEUPSUMMARIES_NORMAL  } from '../../../modules/nf-core/gatk4/gatherpileupsummaries/main'
include { GATK4_GATHERPILEUPSUMMARIES     as GATHERPILEUPSUMMARIES_TUMOR   } from '../../../modules/nf-core/gatk4/gatherpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES        as GETPILEUPSUMMARIES_NORMAL     } from '../../../modules/nf-core/gatk4/getpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES        as GETPILEUPSUMMARIES_TUMOR      } from '../../../modules/nf-core/gatk4/getpileupsummaries/main'
include { GATK4_CALCULATECONTAMINATION    as CALCULATECONTAMINATION        } from '../../../modules/nf-core/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS         as FILTERMUTECTCALLS             } from '../../../modules/nf-core/gatk4/filtermutectcalls/main'


workflow BAM_VARIANT_CALLING_SOMATIC_MUTECT2_MULTI_SAMPLE {
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

    main:
    ch_versions = Channel.empty()
    germline_resource_pileup     = germline_resource_tbi ? germline_resource : Channel.empty()
    germline_resource_pileup_tbi = germline_resource_tbi ?: Channel.empty()


    // Separate normal cram files and remove duplicates
    ch_normal_cram = input.map{ meta, n_cram, n_crai, t_cram, t_crai -> [ meta - meta.subMap('tumor_id') + [id:meta.patient], n_cram, n_crai ] }.unique()
    // Extract tumor cram files
    ch_tumor_cram = input.map{ meta, n_cram, n_crai, t_cram, t_crai -> [ meta - meta.subMap('tumor_id') + [id:meta.patient], t_cram, t_crai ] }
    // Merge normal and tumor crams by patient
    ch_tn_cram = ch_normal_cram.mix(ch_tumor_cram).groupTuple()
    // Combine input and intervals for scatter and gather strategy
    ch_tn_intervals = ch_tn_cram.combine(intervals)
        // Move num_intervals to meta map and reorganize channel for MUTECT2_PAIRED module
        .map{ meta, cram, crai, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], cram, crai, intervals ] 
    }
 
    MUTECT2_PAIRED(
            ch_tn_intervals,
            fasta,
            fai,
            dict,
            germline_resource,
            germline_resource_tbi,
            panel_of_normals,
            panel_of_normals_tbi)


    // Figure out if using intervals or no_intervals
    MUTECT2_PAIRED.out.vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{ mutect2_vcf_branch }

    MUTECT2_PAIRED.out.tbi.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{ mutect2_tbi_branch }

    MUTECT2_PAIRED.out.stats.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{ mutect2_stats_branch }

    MUTECT2_PAIRED.out.f1r2.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{ mutect2_f1r2_branch }

    //Only when using intervals
    MERGE_MUTECT2(
        mutect2_vcf_branch.intervals
        .map{ meta, vcf ->

            new_meta = [id:meta.patient,
                        normal_id:meta.normal_id,
                        num_intervals:meta.num_intervals,
                        patient:meta.patient,
                        sex:meta.sex
                    ]

            [groupKey(new_meta, meta.num_intervals), vcf]
        }.groupTuple(),
        dict
    )

    mutect2_vcf = Channel.empty().mix(
        MERGE_MUTECT2.out.vcf,
        mutect2_vcf_branch.no_intervals)

    mutect2_tbi = Channel.empty().mix(
        MERGE_MUTECT2.out.tbi,
        mutect2_tbi_branch.no_intervals)

    //Merge Mutect2 Stats
    MERGEMUTECTSTATS(
        mutect2_stats_branch.intervals
        .map{ meta, stats ->
            new_meta = [id:meta.patient,
                        normal_id:      meta.normal_id,
                        num_intervals:  meta.num_intervals,
                        patient:        meta.patient,
                        sex:            meta.sex
                    ]

            [groupKey(new_meta, meta.num_intervals), stats]
        }.groupTuple())

    mutect2_stats = Channel.empty().mix(
        MERGEMUTECTSTATS.out.stats,
        mutect2_stats_branch.no_intervals)

    //
    //Generate artifactpriors using learnreadorientationmodel on the f1r2 output of mutect2.
    //
    LEARNREADORIENTATIONMODEL(Channel.empty().mix(
        mutect2_f1r2_branch.intervals
            .map{ meta, f1r2 ->

                new_meta = [
                            id:             meta.patient,
                            normal_id:      meta.normal_id,
                            num_intervals:  meta.num_intervals,
                            patient:        meta.patient,
                            sex:            meta.sex
                        ]

                [groupKey(new_meta, meta.num_intervals), f1r2]
            }.groupTuple(),
        mutect2_f1r2_branch.no_intervals)
    )

    //
    // Generate pileup summary tables using getepileupsummaries. tumor sample should always be passed in as the first input and input list entries of ch_mutect2_in,
    // to ensure correct file order for calculatecontamination.
    input_intervals = input.combine(intervals)
    // Move num_intervals to meta map and reorganize channel for MUTECT2_PAIRED module
        .map{ meta, n_cram, n_crai, t_cram, t_crai, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], n_cram, n_crai, t_cram, t_crai, intervals ] }

    cram_pair_for_pileup = input_intervals.map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals ->
                        [meta, [normal_cram, tumor_cram], [normal_crai, tumor_crai], intervals]
                    }
    pileup = cram_pair_for_pileup.multiMap{  meta, input_list, input_index_list, intervals ->
        tumor: [ meta, input_list[1], input_index_list[1], intervals ]
        normal: [ meta, input_list[0], input_index_list[0], intervals ]
    }

    // Prepare input channel for tumor pileup summaries
    ch_pileup_in_tumor = pileup.tumor.map{meta, cram, crai, intervals ->
                                            [[
                                                id:             meta.tumor_id,
                                                num_intervals:  meta.num_intervals,
                                                patient:        meta.patient,
                                                sex:            meta.sex,
                                                normal_id:      meta.normal_id
                                            ],
                                                cram, crai, intervals]
                                        }.unique()

    GETPILEUPSUMMARIES_TUMOR ( ch_pileup_in_tumor, fasta, fai, dict, germline_resource_pileup, germline_resource_pileup_tbi )

    GETPILEUPSUMMARIES_TUMOR.out.table.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }set{ pileup_table_tumor }

    GATHERPILEUPSUMMARIES_TUMOR(
        pileup_table_tumor.intervals
        .map{ meta, table ->
            new_meta = [
                            id:             meta.id,
                            num_intervals:  meta.num_intervals,
                            patient:        meta.patient,
                            sex:            meta.sex,
                            normal_id:      meta.normal_id
                        ]

            [groupKey(new_meta, meta.num_intervals), table]
        }.groupTuple(),
        dict.map{ meta, dict -> [ dict ] })

    gather_table_tumor = Channel.empty().mix(
        GATHERPILEUPSUMMARIES_TUMOR.out.table,
        pileup_table_tumor.no_intervals).map{ meta, table ->
            new_meta = [
                        id:             meta.id,
                        num_intervals:  meta.num_intervals,
                        patient:        meta.patient,
                        sex:            meta.sex,
                        normal_id:      meta.normal_id
                    ]

            [new_meta, table]
        }

    // Prepare input channel for normal pileup summaries.
    // Remember, the input channel contains tumor-normal pairs, so there will be multiple copies of the normal sample for each tumor for a given patient.
    // Therefore, we use unique function to generate normal pileup summaries once for each patient for better efficiency.
    ch_pileup_in_normal = pileup.normal.map{
                                    meta, cram, crai, intervals ->

                                    [[
                                        id:             meta.normal_id,
                                        num_intervals:  meta.num_intervals,
                                        patient:        meta.patient,
                                        sex:            meta.sex,
                                        normal_id:      meta.normal_id
                                    ],
                                        cram, crai, intervals]
                                }.unique()

    GETPILEUPSUMMARIES_NORMAL ( ch_pileup_in_normal, fasta, fai, dict, germline_resource_pileup, germline_resource_pileup_tbi )

    GETPILEUPSUMMARIES_NORMAL.out.table.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }set{ pileup_table_normal }

    GATHERPILEUPSUMMARIES_NORMAL(
        pileup_table_normal.intervals
        .map{ meta, table ->

            new_meta = [
                            id:             meta.id,
                            num_intervals:  meta.num_intervals,
                            patient:        meta.patient,
                            sex:            meta.sex,
                            normal_id:      meta.normal_id
                        ]

            [groupKey(new_meta, meta.num_intervals), table]
        }.groupTuple(),
        dict.map{ meta, dict -> [ dict ]})

    gather_table_normal = Channel.empty().mix(
        GATHERPILEUPSUMMARIES_NORMAL.out.table,
        pileup_table_normal.no_intervals).map{ meta, table ->

            new_meta = [
                            id:             meta.id,
                            num_intervals:  meta.num_intervals,
                            patient:        meta.patient,
                            sex:            meta.sex,
                            normal_id:      meta.normal_id
                        ]
            [new_meta, table]
        }

    //
    // Contamination and segmentation tables created using calculatecontamination on the pileup summary table.
    // Do some channel magic to generate tumor-normal pairs again.
    // This is necessary because we generated one normal pileup summary for each patient but we need run calculate contamination for each tumor-normal pair.
    ch_tumor_pileup_tables = gather_table_tumor.map{meta, table -> [[id: meta.patient, patient: meta.patient, normal_id: meta.normal_id, sex: meta.sex, num_intervals: meta.num_intervals], meta.id, table]}
    ch_normal_pileup_tables = gather_table_normal.map{meta, table -> [[id: meta.patient, patient: meta.patient, normal_id: meta.normal_id, sex: meta.sex, num_intervals: meta.num_intervals], meta.id, table]}
    ch_calculatecontamination_in_tables = ch_tumor_pileup_tables.combine(ch_normal_pileup_tables, by:0).map{meta, tumor_id, tumor_table, normal_id, normal_table ->
                                                                                                            new_meta = [ id: tumor_id + "_vs_" + normal_id, patient: meta.patient, normal_id: meta.normal_id, sex: meta.sex, num_intervals: meta.num_intervals]
                                                                                                            [new_meta, tumor_table, normal_table]}

    CALCULATECONTAMINATION( ch_calculatecontamination_in_tables)

    //
    // Mutect2 calls filtered by filtermutectcalls using the artifactpriors, contamination and segmentation tables.
    // Reduce the meta to only patient name
    ch_m2_vcf_to_filtermutectcalls = mutect2_vcf.map{meta, vcf -> [[id: meta.patient], vcf]}
    ch_m2_tbi_to_filtermutectcalls = mutect2_tbi.map{meta, tbi -> [[id: meta.patient], tbi]}
    ch_m2_stats_to_filtermutectcalls = mutect2_stats.map{meta, stats -> [[id: meta.patient], stats]}
    ch_artifactpriors_to_filtermutectcalls = LEARNREADORIENTATIONMODEL.out.artifactprior.map{meta, priors -> [[id: meta.patient], priors]}
    ch_seg_to_filtermutectcalls = CALCULATECONTAMINATION.out.segmentation.map{ meta, seg -> [[id: meta.patient], seg]}.groupTuple()
    ch_cont_to_filtermutectcalls = CALCULATECONTAMINATION.out.contamination.map{ meta, cont -> [[id: meta.patient], cont]}.groupTuple()
    // Join the channels by patient name
    ch_filtermutect = ch_m2_vcf_to_filtermutectcalls.join(ch_m2_tbi_to_filtermutectcalls, failOnMismatch: true)
                                                    .join(ch_m2_stats_to_filtermutectcalls, failOnMismatch: true)
                                                    .join(ch_artifactpriors_to_filtermutectcalls, failOnMismatch: true)
                                                    .join(ch_seg_to_filtermutectcalls, failOnMismatch: true)
                                                    .join(ch_cont_to_filtermutectcalls, failOnMismatch: true)

    ch_filtermutect_in = ch_filtermutect.map{ meta, vcf, tbi, stats, orientation, seg, cont -> [meta, vcf, tbi, stats, orientation, seg, cont, []] }

    FILTERMUTECTCALLS ( ch_filtermutect_in, fasta, fai, dict)

    mutect2_vcf_filtered = FILTERMUTECTCALLS.out.vcf.map{ meta, vcf -> [[id:meta.id, variantcaller:"mutect2"], vcf]}
    mutect2_vcf_filtered_tbi = FILTERMUTECTCALLS.out.tbi.map{ meta, tbi -> [[id:meta.id, variantcaller:"mutect2"], tbi]}
    mutect2_vcf_filtered_stats = FILTERMUTECTCALLS.out.stats.map{ meta, stats -> [[id:meta.id, variantcaller:"mutect2"], stats]}

    ch_versions = ch_versions.mix(MUTECT2_PAIRED.out.versions)
    ch_versions = ch_versions.mix(MERGE_MUTECT2.out.versions)
    ch_versions = ch_versions.mix(MERGEMUTECTSTATS.out.versions)
    ch_versions = ch_versions.mix(LEARNREADORIENTATIONMODEL.out.versions)
    ch_versions = ch_versions.mix(GETPILEUPSUMMARIES_NORMAL.out.versions)
    ch_versions = ch_versions.mix(GETPILEUPSUMMARIES_TUMOR.out.versions)
    ch_versions = ch_versions.mix(GATHERPILEUPSUMMARIES_NORMAL.out.versions)
    ch_versions = ch_versions.mix(GATHERPILEUPSUMMARIES_TUMOR.out.versions)
    ch_versions = ch_versions.mix(CALCULATECONTAMINATION.out.versions)
    ch_versions = ch_versions.mix(FILTERMUTECTCALLS.out.versions)

    emit:
    mutect2_vcf            = mutect2_vcf                                    // channel: [ val(meta), [ vcf ] ]
    mutect2_stats          = mutect2_stats                                  // channel: [ val(meta), [ stats ] ]

    artifact_priors        = LEARNREADORIENTATIONMODEL.out.artifactprior    // channel: [ val(meta), [ artifactprior ] ]

    pileup_table_tumor     = gather_table_tumor                             // channel: [ val(meta), [ table_tumor ] ]
    pileup_table_normal    = gather_table_normal                            // channel: [ val(meta), [ table_normal ] ]

    contamination_table    = CALCULATECONTAMINATION.out.contamination       // channel: [ val(meta), [ contamination ] ]
    segmentation_table     = CALCULATECONTAMINATION.out.segmentation        // channel: [ val(meta), [ segmentation ] ]

    filtered_vcf           = mutect2_vcf_filtered                           // channel: [ val(meta), [ vcf ] ]
    filtered_tbi           = mutect2_vcf_filtered_tbi                       // channel: [ val(meta), [ tbi ] ]
    filtered_stats         = mutect2_vcf_filtered_stats                     // channel: [ val(meta), [ stats ] ]

    versions               = ch_versions                                    // channel: [ versions.yml ]
}
