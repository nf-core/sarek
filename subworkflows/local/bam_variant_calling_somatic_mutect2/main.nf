//
// Run GATK mutect2 in tumor normal mode, getepileupsummaries, calculatecontamination, learnreadorientationmodel and filtermutectcalls
//

include { GATK4_MERGEVCFS                 as MERGE_MUTECT2               } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_CALCULATECONTAMINATION    as CALCULATECONTAMINATION      } from '../../../modules/nf-core/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS         as FILTERMUTECTCALLS           } from '../../../modules/nf-core/gatk4/filtermutectcalls/main'
include { GATK4_GATHERPILEUPSUMMARIES     as GATHERPILEUPSUMMARIES_NORMAL} from '../../../modules/nf-core/gatk4/gatherpileupsummaries/main'
include { GATK4_GATHERPILEUPSUMMARIES     as GATHERPILEUPSUMMARIES_TUMOR } from '../../../modules/nf-core/gatk4/gatherpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES        as GETPILEUPSUMMARIES_NORMAL   } from '../../../modules/nf-core/gatk4/getpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES        as GETPILEUPSUMMARIES_TUMOR    } from '../../../modules/nf-core/gatk4/getpileupsummaries/main'
include { GATK4_LEARNREADORIENTATIONMODEL as LEARNREADORIENTATIONMODEL   } from '../../../modules/nf-core/gatk4/learnreadorientationmodel/main'
include { GATK4_MERGEMUTECTSTATS          as MERGEMUTECTSTATS            } from '../../../modules/nf-core/gatk4/mergemutectstats/main'
include { GATK4_MUTECT2                   as MUTECT2_PAIRED              } from '../../../modules/nf-core/gatk4/mutect2/main'

workflow BAM_VARIANT_CALLING_SOMATIC_MUTECT2 {
    take:
    input                     // channel: [ val(meta), [ input ], [ input_index ], [which_norm] ]
    fasta                     // channel: /path/to/reference/fasta
    fai                       // channel: /path/to/reference/fasta/index
    dict                      // channel: /path/to/reference/fasta/dictionary
    germline_resource         // channel: /path/to/germline/resource
    germline_resource_tbi     // channel: /path/to/germline/index
    panel_of_normals          // channel: /path/to/panel/of/normals
    panel_of_normals_tbi      // channel: /path/to/panel/of/normals/index

    main:
    ch_versions = Channel.empty()

    //
    //Perform variant calling using mutect2 module in tumor single mode.
    //
    MUTECT2_PAIRED(input,
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

            new_meta = [
                        id:meta.tumor_id + "_vs_" + meta.normal_id,
                        normal_id:meta.normal_id,
                        num_intervals:meta.num_intervals,
                        patient:meta.patient,
                        sex:meta.sex,
                        tumor_id:meta.tumor_id
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

            new_meta = [
                        id:             meta.tumor_id + "_vs_" + meta.normal_id,
                        normal_id:      meta.normal_id,
                        num_intervals:  meta.num_intervals,
                        patient:        meta.patient,
                        sex:            meta.sex,
                        tumor_id:       meta.tumor_id
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
                            id:             meta.tumor_id + "_vs_" + meta.normal_id,
                            normal_id:      meta.normal_id,
                            num_intervals:  meta.num_intervals,
                            patient:        meta.patient,
                            sex:            meta.sex,
                            tumor_id:       meta.tumor_id,
                        ]

                [groupKey(new_meta, meta.num_intervals), f1r2]
            }.groupTuple(),
        mutect2_f1r2_branch.no_intervals)
    )

    //
    //Generate pileup summary tables using getepileupsummaries. tumor sample should always be passed in as the first input and input list entries of ch_mutect2_in,
    //to ensure correct file order for calculatecontamination.
    pileup = input.multiMap{  meta, input_list, input_index_list, intervals ->
        tumor: [ meta, input_list[1], input_index_list[1], intervals ]
        normal: [ meta, input_list[0], input_index_list[0], intervals ]
    }

    germline_resource_pileup = germline_resource_tbi ? germline_resource : Channel.empty()
    germline_resource_pileup_tbi = germline_resource_tbi ?: Channel.empty()
    GETPILEUPSUMMARIES_TUMOR ( pileup.tumor.map{
                                    meta, cram, crai, intervals ->

                                    [[
                                        id:             meta.tumor_id,
                                        normal_id:      meta.normal_id,
                                        num_intervals:  meta.num_intervals,
                                        patient:        meta.patient,
                                        sex:            meta.sex,
                                        tumor_id:       meta.tumor_id,
                                    ],
                                        cram, crai, intervals]
                                },
                                fasta, fai, dict, germline_resource_pileup, germline_resource_pileup_tbi )

    GETPILEUPSUMMARIES_NORMAL ( pileup.normal.map{
                                    meta, cram, crai, intervals ->

                                    [[
                                        id:             meta.normal_id,
                                        normal_id:      meta.normal_id,
                                        num_intervals:  meta.num_intervals,
                                        patient:        meta.patient,
                                        sex:            meta.sex,
                                        tumor_id:       meta.tumor_id,
                                    ],
                                        cram, crai, intervals]
                                },
                                fasta, fai, dict, germline_resource_pileup, germline_resource_pileup_tbi )

    GETPILEUPSUMMARIES_NORMAL.out.table.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }set{ pileup_table_normal }

    GETPILEUPSUMMARIES_TUMOR.out.table.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }set{ pileup_table_tumor }

    //Merge Pileup Summaries
    GATHERPILEUPSUMMARIES_NORMAL(
        GETPILEUPSUMMARIES_NORMAL.out.table
        .map{ meta, table ->

            new_meta = [
                            id:             meta.normal_id,
                            normal_id:      meta.normal_id,
                            num_intervals:  meta.num_intervals,
                            patient:        meta.patient,
                            sex:            meta.sex,
                            tumor_id:       meta.tumor_id,
                        ]

            [groupKey(new_meta, meta.num_intervals), table]
        }.groupTuple(),
        dict)

    gather_table_normal = Channel.empty().mix(
        GATHERPILEUPSUMMARIES_NORMAL.out.table,
        pileup_table_normal.no_intervals).map{ meta, table ->

            new_meta = [
                            id:             meta.tumor_id + "_vs_" + meta.normal_id,
                            normal_id:      meta.normal_id,
                            num_intervals:  meta.num_intervals,
                            patient:        meta.patient,
                            sex:            meta.sex,
                            tumor_id:       meta.tumor_id,
                        ]
            [new_meta, table]
        }

    GATHERPILEUPSUMMARIES_TUMOR(
        GETPILEUPSUMMARIES_TUMOR.out.table
        .map{ meta, table ->
            new_meta = [
                            id:             meta.tumor_id,
                            normal_id:      meta.normal_id,
                            num_intervals:  meta.num_intervals,
                            patient:        meta.patient,
                            sex:            meta.sex,
                            tumor_id:       meta.tumor_id,
                        ]

            [groupKey(new_meta, meta.num_intervals), table]
        }.groupTuple(),
        dict)

    gather_table_tumor = Channel.empty().mix(
        GATHERPILEUPSUMMARIES_TUMOR.out.table,
        pileup_table_tumor.no_intervals).map{ meta, table ->
            new_meta = [
                        id:             meta.tumor_id + "_vs_" + meta.normal_id,
                        normal_id:      meta.normal_id,
                        num_intervals:  meta.num_intervals,
                        patient:        meta.patient,
                        sex:            meta.sex,
                        tumor_id:       meta.tumor_id,
                    ]

            [new_meta, table]
        }

    //
    //Contamination and segmentation tables created using calculatecontamination on the pileup summary table.
    //
    CALCULATECONTAMINATION ( gather_table_tumor.join(gather_table_normal) )

    //
    //Mutect2 calls filtered by filtermutectcalls using the artifactpriors, contamination and segmentation tables.
    //
    ch_filtermutect    = mutect2_vcf.join(mutect2_tbi)
                                    .join(mutect2_stats)
                                    .join(LEARNREADORIENTATIONMODEL.out.artifactprior)
                                    .join(CALCULATECONTAMINATION.out.segmentation)
                                    .join(CALCULATECONTAMINATION.out.contamination)
    ch_filtermutect_in = ch_filtermutect.map{ meta, vcf, tbi, stats, orientation, seg, cont -> [meta, vcf, tbi, stats, orientation, seg, cont, []] }

    FILTERMUTECTCALLS ( ch_filtermutect_in, fasta, fai, dict )

    ch_versions = ch_versions.mix(MERGE_MUTECT2.out.versions)
    ch_versions = ch_versions.mix(CALCULATECONTAMINATION.out.versions)
    ch_versions = ch_versions.mix(FILTERMUTECTCALLS.out.versions)
    ch_versions = ch_versions.mix(GETPILEUPSUMMARIES_NORMAL.out.versions)
    ch_versions = ch_versions.mix(GETPILEUPSUMMARIES_TUMOR.out.versions)
    ch_versions = ch_versions.mix(GATHERPILEUPSUMMARIES_NORMAL.out.versions)
    ch_versions = ch_versions.mix(GATHERPILEUPSUMMARIES_TUMOR.out.versions)
    ch_versions = ch_versions.mix(LEARNREADORIENTATIONMODEL.out.versions)
    ch_versions = ch_versions.mix(MERGEMUTECTSTATS.out.versions)
    ch_versions = ch_versions.mix(MUTECT2_PAIRED.out.versions)

    emit:
    mutect2_vcf            = mutect2_vcf                                    // channel: [ val(meta), [ vcf ] ]
    mutect2_stats          = mutect2_stats                                  // channel: [ val(meta), [ stats ] ]

    artifact_priors        = LEARNREADORIENTATIONMODEL.out.artifactprior    // channel: [ val(meta), [ artifactprior ] ]

    pileup_table_tumor     = gather_table_tumor                             // channel: [ val(meta), [ table_tumor ] ]
    pileup_table_normal    = gather_table_normal                            // channel: [ val(meta), [ table_normal ] ]

    contamination_table    = CALCULATECONTAMINATION.out.contamination       // channel: [ val(meta), [ contamination ] ]
    segmentation_table     = CALCULATECONTAMINATION.out.segmentation        // channel: [ val(meta), [ segmentation ] ]

    filtered_vcf           = FILTERMUTECTCALLS.out.vcf.map{ meta, vcf -> [[patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, sex:meta.sex, id:meta.tumor_id + "_vs_" + meta.normal_id, num_intervals:meta.num_intervals, variantcaller:"mutect2"],
                                                                            vcf]} // channel: [ val(meta), [ vcf ] ]
    filtered_tbi           = FILTERMUTECTCALLS.out.tbi                      // channel: [ val(meta), [ tbi ] ]
    filtered_stats         = FILTERMUTECTCALLS.out.stats                    // channel: [ val(meta), [ stats ] ]

    versions               = ch_versions                                    // channel: [ versions.yml ]
}
