//
// Run GATK mutect2 in tumor normal mode, getepileupsummaries, calculatecontamination, learnreadorientationmodel and filtermutectcalls
//

include { BGZIP                           as BGZIP_MUTECT2               } from '../../../modules/local/bgzip'
include { CONCAT_VCF                      as CONCAT_VCF_MUTECT2          } from '../../../modules/local/concat_vcf/main'

include { GATK4_MUTECT2                   as MUTECT2                     } from '../../../modules/nf-core/modules/gatk4/mutect2/main'
include { GATK4_MERGEMUTECTSTATS          as MERGEMUTECTSTATS            } from '../../../modules/local/gatk4/mergemutectstats'
include { GATK4_LEARNREADORIENTATIONMODEL as LEARNREADORIENTATIONMODEL   } from '../../../modules/nf-core/modules/gatk4/learnreadorientationmodel/main'
include { GATK4_GATHERPILEUPSUMMARIES     as GATHERPILEUPSUMMARIES_TUMOR } from '../../../modules/local/gatk4/gatherpileupsummaries'
include { GATK4_GATHERPILEUPSUMMARIES     as GATHERPILEUPSUMMARIES_NORMAL} from '../../../modules/local/gatk4/gatherpileupsummaries'
include { GATK4_GETPILEUPSUMMARIES        as GETPILEUPSUMMARIES_TUMOR    } from '../../../modules/nf-core/modules/gatk4/getpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES        as GETPILEUPSUMMARIES_NORMAL   } from '../../../modules/nf-core/modules/gatk4/getpileupsummaries/main'
include { GATK4_CALCULATECONTAMINATION    as CALCULATECONTAMINATION      } from '../../../modules/nf-core/modules/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS         as FILTERMUTECTCALLS           } from '../../../modules/nf-core/modules/gatk4/filtermutectcalls/main'

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
    no_intervals
    num_intervals
    intervals_bed_combine_gz

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

    pileup_tumor_input = input.map {
        meta, input_list, input_index_list, intervals, which_norm ->
        tumor_id = meta.tumor_id
        id = intervals ? tumor_id + "_" + intervals.baseName : tumor_id
        [[id: id, normal_id: meta.normal_id, tumor_id: meta.tumor_id, gender: meta.gender, patient: meta.patient ], input_list[1], input_index_list[1], intervals]
    }

    pileup_normal_input = input.map {
        meta, input_list, input_index_list, intervals, which_norm ->
        normal_id = meta.normal_id
        id = intervals ? normal_id + "_" + intervals.baseName : normal_id
        [[id: id, normal_id: meta.normal_id, tumor_id: meta.tumor_id, gender: meta.gender, patient: meta.patient ], input_list[0], input_index_list[0], intervals]
    }
    GETPILEUPSUMMARIES_TUMOR ( pileup_tumor_input, fasta, fai, dict, germline_resource, germline_resource_tbi )
    GETPILEUPSUMMARIES_NORMAL ( pileup_normal_input, fasta, fai, dict, germline_resource, germline_resource_tbi )
    ch_versions = ch_versions.mix(GETPILEUPSUMMARIES_NORMAL.out.versions)

    if(no_intervals){
        mutect2_vcf_gz_tbi = MUTECT2.out.vcf.join(MUTECT2.out.tbi)
        mutect2_stats = MUTECT2.out.stats
        pileup_table_tumor= GETPILEUPSUMMARIES_TUMOR.out.table
        pileup_table_normal= GETPILEUPSUMMARIES_NORMAL.out.table

    }else{

        //Merge Mutect2 VCF
        BGZIP_MUTECT2(MUTECT2.out.vcf)

        BGZIP_MUTECT2.out.vcf.map{ meta, vcf ->
            id = meta.tumor_id + "_vs_" + meta.normal_id
            [[id: id, normal_id: meta.normal_id, tumor_id: meta.tumor_id, gender: meta.gender, patient: meta.patient ], vcf]
        }.set{bgzip_mutect2}

        mutect2_vcf_to_concat = bgzip_mutect2.groupTuple(size: num_intervals)

        CONCAT_VCF_MUTECT2(mutect2_vcf_to_concat, fai, intervals_bed_combine_gz)
        mutect2_vcf_gz_tbi = CONCAT_VCF_MUTECT2.out.vcf.join(CONCAT_VCF_MUTECT2.out.tbi)

        ch_versions = ch_versions.mix(BGZIP_MUTECT2.out.versions)
        ch_versions = ch_versions.mix(CONCAT_VCF_MUTECT2.out.versions)

        //Merge Muteect2 Stats
        MUTECT2.out.stats.map{ meta, stats ->
            id = meta.tumor_id + "_vs_" + meta.normal_id
            [[id: id, normal_id: meta.normal_id, tumor_id: meta.tumor_id, gender: meta.gender, patient: meta.patient ], stats]
        }.groupTuple(size: num_intervals).set{mutect2_stats_to_merge}

        MERGEMUTECTSTATS(mutect2_stats_to_merge)
        mutect2_stats = MERGEMUTECTSTATS.out.stats
        ch_versions = ch_versions.mix(MERGEMUTECTSTATS.out.versions)

        //Merge Pileup Summaries
        pileup_tumor_tables_to_gather = GETPILEUPSUMMARIES_TUMOR.out.table.map{ meta, table ->
            [[id: meta.tumor_id, normal_id: meta.normal_id, tumor_id: meta.tumor_id, gender: meta.gender, patient: meta.patient ], table]
        }.groupTuple(size: num_intervals)

        GATHERPILEUPSUMMARIES_TUMOR(pileup_tumor_tables_to_gather, dict)
        GATHERPILEUPSUMMARIES_TUMOR.out.table.map{ meta, table ->
            id = meta.tumor_id + "_vs_" + meta.normal_id
            [[id: id, normal_id: meta.normal_id, tumor_id: meta.tumor_id, gender: meta.gender, patient: meta.patient ], table]
        }.set{pileup_table_tumor}

        pileup_normal_tables_to_gather = GETPILEUPSUMMARIES_NORMAL.out.table.map{ meta, table ->
            [[id: meta.normal_id, normal_id: meta.normal_id, tumor_id: meta.tumor_id, gender: meta.gender, patient: meta.patient ], table]
        }.groupTuple(size: num_intervals)

        GATHERPILEUPSUMMARIES_NORMAL(pileup_normal_tables_to_gather, dict)
        GATHERPILEUPSUMMARIES_NORMAL.out.table.map{ meta, table ->
            id = meta.tumor_id + "_vs_" + meta.normal_id
            [[id: id, normal_id: meta.normal_id, tumor_id: meta.tumor_id, gender: meta.gender, patient: meta.patient ], table]
        }.set{pileup_table_normal}

    }

    //
    //Generate artifactpriors using learnreadorientationmodel on the f1r2 output of mutect2.
    //
    MUTECT2.out.f1r2.map{ meta, f1f2 ->
        id = meta.tumor_id + "_vs_" + meta.normal_id
        [[id: id, normal_id: meta.normal_id, tumor_id: meta.tumor_id, gender: meta.gender, patient: meta.patient ], f1f2]
    }.groupTuple(size: num_intervals)
    .set{ch_learnread_in}

    LEARNREADORIENTATIONMODEL (ch_learnread_in)

    ch_versions = ch_versions.mix(LEARNREADORIENTATIONMODEL.out.versions)

    //
    //Contamination and segmentation tables created using calculatecontamination on the pileup summary table.
    //
    ch_calccon_in       = pileup_table_tumor.join(pileup_table_normal)
    CALCULATECONTAMINATION ( ch_calccon_in, true )
    ch_versions   = ch_versions.mix(CALCULATECONTAMINATION.out.versions)

    //
    //Mutect2 calls filtered by filtermutectcalls using the artifactpriors, contamination and segmentation tables.
    //
    ch_filtermutect      = mutect2_vcf_gz_tbi.join(mutect2_stats)
                                            .join(LEARNREADORIENTATIONMODEL.out.artifactprior)
                                            .join(CALCULATECONTAMINATION.out.segmentation)
                                            .join(CALCULATECONTAMINATION.out.contamination)
    ch_filtermutect.map{ meta, vcf, tbi, stats, orientation, seg, cont ->
        [meta, vcf, tbi, stats, orientation, seg, cont, []]
    }.set{ch_filtermutect_in}

    FILTERMUTECTCALLS ( ch_filtermutect_in, fasta, fai, dict )
    ch_versions              = ch_versions.mix(FILTERMUTECTCALLS.out.versions)

    emit:
    mutect2_vcf_gz_tbi     = mutect2_vcf_gz_tbi                             // channel: [ val(meta), [ vcf ] ]
    mutect2_stats          = mutect2_stats                              // channel: [ val(meta), [ stats ] ]
    //mutect2_f1r2           = MUTECT2.out.f1r2                               // channel: [ val(meta), [ f1r2 ] ]

    artifact_priors        = LEARNREADORIENTATIONMODEL.out.artifactprior // channel: [ val(meta), [ artifactprior ] ]

    pileup_table_tumor     = pileup_table_tumor          // channel: [ val(meta), [ table_tumor ] ]
    pileup_table_normal    = pileup_table_normal         // channel: [ val(meta), [ table_normal ] ]

    contamination_table    = CALCULATECONTAMINATION.out.contamination    // channel: [ val(meta), [ contamination ] ]
    segmentation_table     = CALCULATECONTAMINATION.out.segmentation     // channel: [ val(meta), [ segmentation ] ]

    filtered_vcf           = FILTERMUTECTCALLS.out.vcf                   // channel: [ val(meta), [ vcf ] ]
    filtered_tbi           = FILTERMUTECTCALLS.out.tbi                   // channel: [ val(meta), [ tbi ] ]
    filtered_stats         = FILTERMUTECTCALLS.out.stats                 // channel: [ val(meta), [ stats ] ]

    versions               = ch_versions                                           // channel: [ versions.yml ]
}
