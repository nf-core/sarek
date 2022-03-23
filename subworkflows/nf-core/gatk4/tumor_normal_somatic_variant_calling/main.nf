//
// Run GATK mutect2 in tumor normal mode, getepileupsummaries, calculatecontamination, learnreadorientationmodel and filtermutectcalls
//

include { BGZIP                           as BGZIP_MUTECT2               } from '../../../../modules/local/bgzip'
include { CONCAT_VCF                      as CONCAT_MUTECT2              } from '../../../../modules/local/concat_vcf/main'
include { GATK4_MUTECT2                   as MUTECT2                     } from '../../../../modules/nf-core/modules/gatk4/mutect2/main'
include { GATK4_MERGEMUTECTSTATS          as MERGEMUTECTSTATS            } from '../../../../modules/nf-core/modules/gatk4/mergemutectstats/main'
include { GATK4_LEARNREADORIENTATIONMODEL as LEARNREADORIENTATIONMODEL   } from '../../../../modules/nf-core/modules/gatk4/learnreadorientationmodel/main'
include { GATK4_GATHERPILEUPSUMMARIES     as GATHERPILEUPSUMMARIES_TUMOR } from '../../../../modules/nf-core/modules/gatk4/gatherpileupsummaries/main'
include { GATK4_GATHERPILEUPSUMMARIES     as GATHERPILEUPSUMMARIES_NORMAL} from '../../../../modules/nf-core/modules/gatk4/gatherpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES        as GETPILEUPSUMMARIES_TUMOR    } from '../../../../modules/nf-core/modules/gatk4/getpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES        as GETPILEUPSUMMARIES_NORMAL   } from '../../../../modules/nf-core/modules/gatk4/getpileupsummaries/main'
include { GATK4_CALCULATECONTAMINATION    as CALCULATECONTAMINATION      } from '../../../../modules/nf-core/modules/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS         as FILTERMUTECTCALLS           } from '../../../../modules/nf-core/modules/gatk4/filtermutectcalls/main'

workflow GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING {
    take:
    input                     // channel: [ val(meta), [ input ], [ input_index ], [which_norm] ]
    fasta                     // channel: /path/to/reference/fasta
    fai                       // channel: /path/to/reference/fasta/index
    dict                      // channel: /path/to/reference/fasta/dictionary
    germline_resource         // channel: /path/to/germline/resource
    germline_resource_tbi     // channel: /path/to/germline/index
    panel_of_normals          // channel: /path/to/panel/of/normals
    panel_of_normals_tbi      // channel: /path/to/panel/of/normals/index
    intervals_bed_combine_gz
    num_intervals

    main:
    ch_versions = Channel.empty()

    //
    //Perform variant calling using mutect2 module in tumor single mode.
    //
    MUTECT2 ( input, false, false, false, fasta, fai, dict, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi )
    ch_versions = ch_versions.mix(MUTECT2.out.versions)

    //
    //Generate pileup summary tables using getepileupsummaries. tumor sample should always be passed in as the first input and input list entries of ch_mutect2_in,
    //to ensure correct file order for calculatecontamination.
    input.multiMap{  meta, input_list, input_index_list, intervals, which_norm ->
        tumor: [ meta, input_list[1], input_index_list[1], intervals ]
        normal: [ meta, input_list[0], input_index_list[0], intervals ]
    }.set{pileup}

    GETPILEUPSUMMARIES_TUMOR ( pileup.tumor, fasta, fai, dict, germline_resource, germline_resource_tbi )
    GETPILEUPSUMMARIES_NORMAL ( pileup.normal, fasta, fai, dict, germline_resource, germline_resource_tbi )
    ch_versions = ch_versions.mix(GETPILEUPSUMMARIES_NORMAL.out.versions)

    // Figure out if using intervals or no_intervals
    MUTECT2.out.vcf.branch{
            intervals:    num_intervals > 1
            no_intervals: num_intervals == 1
        }.set{ mutect2_vcf }

    MUTECT2.out.tbi.branch{
            intervals:    num_intervals > 1
            no_intervals: num_intervals == 1
        }.set{ mutect2_tbi }

    MUTECT2.out.stats.branch{
            intervals:    num_intervals > 1
            no_intervals: num_intervals == 1
        }.set{ mutect2_stats }

    GETPILEUPSUMMARIES_NORMAL.out.table.branch{
            intervals:    num_intervals > 1
            no_intervals: num_intervals == 1
        }set{ pileup_table_normal }

    GETPILEUPSUMMARIES_TUMOR.out.table.branch{
            intervals:    num_intervals > 1
            no_intervals: num_intervals == 1
        }set{ pileup_table_tumor }


    //Only when using intervals

    //Merge Mutect2 VCF
    BGZIP_MUTECT2(MUTECT2.out.vcf)

    CONCAT_MUTECT2(
        BGZIP_MUTECT2.out.vcf.map{ meta, vcf ->
            [[id: meta.tumor_id + "_vs_" + meta.normal_id, normal_id: meta.normal_id, tumor_id: meta.tumor_id, gender: meta.gender, patient: meta.patient ], vcf]
        }.groupTuple(size: num_intervals),
        fai,
        intervals_bed_combine_gz)

    mutect2_vcf = Channel.empty().mix(
        CONCAT_MUTECT2.out.vcf,
        mutect2_vcf.no_intervals)

    mutect2_tbi = Channel.empty().mix(
        CONCAT_MUTECT2.out.tbi,
        mutect2_tbi.no_intervals)

    ch_versions = ch_versions.mix(BGZIP_MUTECT2.out.versions)
    ch_versions = ch_versions.mix(CONCAT_MUTECT2.out.versions)

    //Merge Muteect2 Stats
    MERGEMUTECTSTATS(mutect2_stats.intervals.map{ meta, stats ->
        [[id: meta.tumor_id + "_vs_" + meta.normal_id, normal_id: meta.normal_id, tumor_id: meta.tumor_id, gender: meta.gender, patient: meta.patient ], stats]
    }.groupTuple(size: num_intervals))

    mutect2_stats = Channel.empty().mix(
        MERGEMUTECTSTATS.out.stats,
        mutect2_stats.no_intervals)

    //Merge Pileup Summaries
    GATHERPILEUPSUMMARIES_NORMAL(
        GETPILEUPSUMMARIES_NORMAL.out.table.map{ meta, table ->
            [[id: meta.normal_id, normal_id: meta.normal_id, tumor_id: meta.tumor_id, gender: meta.gender, patient: meta.patient ], table]
        }.groupTuple(size: num_intervals),
        dict)

    gather_table_normal = Channel.empty().mix(
        GATHERPILEUPSUMMARIES_NORMAL.out.table.map{ meta, table ->
            [[id: meta.tumor_id + "_vs_" + meta.normal_id, normal_id: meta.normal_id, tumor_id: meta.tumor_id, gender: meta.gender, patient: meta.patient ], table]
        },
        pileup_table_normal.no_intervals)

    GATHERPILEUPSUMMARIES_TUMOR( GETPILEUPSUMMARIES_TUMOR.out.table.map{ meta, table ->
            [[id: meta.tumor_id, normal_id: meta.normal_id, tumor_id: meta.tumor_id, gender: meta.gender, patient: meta.patient ], table]
        }.groupTuple(size: num_intervals),
        dict)

    gather_table_tumor = Channel.empty().mix(
        GATHERPILEUPSUMMARIES_TUMOR.out.table.map{ meta, table ->
            [[id: meta.tumor_id + "_vs_" + meta.normal_id, normal_id: meta.normal_id, tumor_id: meta.tumor_id, gender: meta.gender, patient: meta.patient ], table]
        },
        pileup_table_tumor.no_intervals)

    //
    //Generate artifactpriors using learnreadorientationmodel on the f1r2 output of mutect2.
    //
    MUTECT2.out.f1r2.map{ meta, f1f2 ->
        [[id: meta.tumor_id + "_vs_" + meta.normal_id, normal_id: meta.normal_id, tumor_id: meta.tumor_id, gender: meta.gender, patient: meta.patient ], f1f2]
    }.groupTuple(size: num_intervals)
    .set{ch_learnread_in}

    LEARNREADORIENTATIONMODEL (ch_learnread_in)
    ch_versions = ch_versions.mix(LEARNREADORIENTATIONMODEL.out.versions)

    //
    //Contamination and segmentation tables created using calculatecontamination on the pileup summary table.
    //
    ch_calccon_in       = gather_table_tumor.join(gather_table_normal)
    CALCULATECONTAMINATION ( ch_calccon_in, true )
    ch_versions   = ch_versions.mix(CALCULATECONTAMINATION.out.versions)

    //
    //Mutect2 calls filtered by filtermutectcalls using the artifactpriors, contamination and segmentation tables.
    //
    ch_filtermutect      = mutect2_vcf.join(mutect2_tbi)
                                      .join(mutect2_stats)
                                      .join(LEARNREADORIENTATIONMODEL.out.artifactprior)
                                      .join(CALCULATECONTAMINATION.out.segmentation)
                                      .join(CALCULATECONTAMINATION.out.contamination)
    ch_filtermutect.map{ meta, vcf, tbi, stats, orientation, seg, cont ->
        [meta, vcf, tbi, stats, orientation, seg, cont, []]
    }.set{ch_filtermutect_in}

    FILTERMUTECTCALLS ( ch_filtermutect_in, fasta, fai, dict )
    ch_versions              = ch_versions.mix(FILTERMUTECTCALLS.out.versions)

    emit:
    mutect2_vcf            = mutect2_vcf                                    // channel: [ val(meta), [ vcf ] ]
    mutect2_stats          = mutect2_stats                                  // channel: [ val(meta), [ stats ] ]
    mutect2_f1r2           = MUTECT2.out.f1r2                               // channel: [ val(meta), [ f1r2 ] ]

    artifact_priors        = LEARNREADORIENTATIONMODEL.out.artifactprior    // channel: [ val(meta), [ artifactprior ] ]

    pileup_table_tumor     = gather_table_tumor                             // channel: [ val(meta), [ table_tumor ] ]
    pileup_table_normal    = gather_table_normal                            // channel: [ val(meta), [ table_normal ] ]

    contamination_table    = CALCULATECONTAMINATION.out.contamination       // channel: [ val(meta), [ contamination ] ]
    segmentation_table     = CALCULATECONTAMINATION.out.segmentation        // channel: [ val(meta), [ segmentation ] ]

    filtered_vcf           = FILTERMUTECTCALLS.out.vcf                      // channel: [ val(meta), [ vcf ] ]
    filtered_tbi           = FILTERMUTECTCALLS.out.tbi                      // channel: [ val(meta), [ tbi ] ]
    filtered_stats         = FILTERMUTECTCALLS.out.stats                    // channel: [ val(meta), [ stats ] ]

    versions               = ch_versions                                    // channel: [ versions.yml ]
}
