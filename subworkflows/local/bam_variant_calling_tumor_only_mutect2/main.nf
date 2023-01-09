//
// Run GATK mutect2 in tumor only mode, getepileupsummaries, calculatecontamination and filtermutectcalls
//

include { GATK4_MERGEVCFS                 as MERGE_MUTECT2             } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_CALCULATECONTAMINATION    as CALCULATECONTAMINATION    } from '../../../modules/nf-core/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS         as FILTERMUTECTCALLS         } from '../../../modules/nf-core/gatk4/filtermutectcalls/main'
include { GATK4_GETPILEUPSUMMARIES        as GETPILEUPSUMMARIES        } from '../../../modules/nf-core/gatk4/getpileupsummaries/main'
include { GATK4_GATHERPILEUPSUMMARIES     as GATHERPILEUPSUMMARIES     } from '../../../modules/nf-core/gatk4/gatherpileupsummaries/main'
include { GATK4_LEARNREADORIENTATIONMODEL as LEARNREADORIENTATIONMODEL } from '../../../modules/nf-core/gatk4/learnreadorientationmodel/main'
include { GATK4_MERGEMUTECTSTATS          as MERGEMUTECTSTATS          } from '../../../modules/nf-core/gatk4/mergemutectstats/main'
include { GATK4_MUTECT2                   as MUTECT2                   } from '../../../modules/nf-core/gatk4/mutect2/main'

workflow BAM_VARIANT_CALLING_TUMOR_ONLY_MUTECT2 {
    take:
    input                     // channel: [ val(meta), [ input ], [ input_index ], [intervals], [] ]
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
    MUTECT2(input,
            fasta,
            fai,
            dict,
            germline_resource,
            germline_resource_tbi,
            panel_of_normals,
            panel_of_normals_tbi)

    // Figure out if using intervals or no_intervals
    MUTECT2.out.vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{ vcf_branch }

    MUTECT2.out.tbi.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{ mutect2_tbi_branch }

    MUTECT2.out.stats.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{ stats_branch }

    MUTECT2.out.f1r2.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{ mutect2_f1r2_branch }

    //Only when using intervals
    //Merge Mutect2 VCF

    MERGE_MUTECT2(
        vcf_branch.intervals
        .map{ meta, vcf ->
            new_meta = [
                        id:             meta.sample,
                        num_intervals:  meta.num_intervals,
                        patient:        meta.patient,
                        sample:         meta.sample,
                        sex:            meta.sex,
                        status:         meta.status,
                    ]

            [groupKey(new_meta, meta.num_intervals), vcf]
        }.groupTuple(),
        dict.map{ it -> [ [ id:'dict' ], it ] })

    vcf = Channel.empty().mix(
        MERGE_MUTECT2.out.vcf,
        vcf_branch.no_intervals)

    mutect2_tbi = Channel.empty().mix(
        MERGE_MUTECT2.out.tbi,
        mutect2_tbi_branch.no_intervals)

    //Merge Mutect2 Stats
    MERGEMUTECTSTATS(
        stats_branch.intervals
        .map{ meta, stats ->
            new_meta = [
                        id:             meta.sample,
                        num_intervals:  meta.num_intervals,
                        patient:        meta.patient,
                        sample:         meta.sample,
                        sex:            meta.sex,
                        status:         meta.status,
                    ]

            [groupKey(new_meta, meta.num_intervals), stats]
        }.groupTuple())

    stats = Channel.empty().mix(
        MERGEMUTECTSTATS.out.stats,
        stats_branch.no_intervals)

    //
    //Generate artifactpriors using learnreadorientationmodel on the f1r2 output of mutect2.
    //
    LEARNREADORIENTATIONMODEL(
        Channel.empty().mix(
            mutect2_f1r2_branch.intervals
            .map{ meta, f1r2 ->
                new_meta = [
                            id:             meta.sample,
                            num_intervals:  meta.num_intervals,
                            patient:        meta.patient,
                            sample:         meta.sample,
                            sex:            meta.sex,
                            status:         meta.status,
                        ]

                [groupKey(new_meta, meta.num_intervals), f1r2]
            }.groupTuple(),
            mutect2_f1r2_branch.no_intervals))

    //
    //Generate pileup summary table using getepileupsummaries.
    //
    germline_resource_pileup = germline_resource_tbi ? germline_resource : Channel.empty()
    germline_resource_pileup_tbi = germline_resource_tbi ?: Channel.empty()
    GETPILEUPSUMMARIES(input , fasta, fai, dict, germline_resource_pileup , germline_resource_pileup_tbi)

    GETPILEUPSUMMARIES.out.table.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }set{ pileup_table_branch }

    //Merge Pileup Summaries
    GATHERPILEUPSUMMARIES(
        GETPILEUPSUMMARIES.out.table
        .map{ meta, table ->
            new_meta = [
                        id:             meta.sample,
                        num_intervals:  meta.num_intervals,
                        patient:        meta.patient,
                        sample:         meta.sample,
                        sex:            meta.sex,
                        status:         meta.status,
                    ]

            [groupKey(new_meta, meta.num_intervals), table]
        }.groupTuple(),
        dict)

    pileup_table = Channel.empty().mix(
        GATHERPILEUPSUMMARIES.out.table,
        pileup_table_branch.no_intervals)

    //
    //Contamination and segmentation tables created using calculatecontamination on the pileup summary table.
    //
    table_contamination = pileup_table.map{ meta, table -> [ meta, table, [] ] }
    CALCULATECONTAMINATION(table_contamination)

    //
    //Mutect2 calls filtered by filtermutectcalls using the contamination and segmentation tables.
    //
    filtermutect = vcf.join(mutect2_tbi)
        .join(stats)
        .join(LEARNREADORIENTATIONMODEL.out.artifactprior)
        .join(CALCULATECONTAMINATION.out.segmentation)
        .join(CALCULATECONTAMINATION.out.contamination)
        .map{ meta, vcf, tbi, stats, artifactprior, seg, cont -> [ meta, vcf, tbi, stats, artifactprior, seg, cont, [] ] }

    FILTERMUTECTCALLS(filtermutect, fasta, fai, dict)

    vcf_filtered = FILTERMUTECTCALLS.out.vcf
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'mutect2' ], vcf ] }

    ch_versions = ch_versions.mix(MERGE_MUTECT2.out.versions)
    ch_versions = ch_versions.mix(CALCULATECONTAMINATION.out.versions)
    ch_versions = ch_versions.mix(FILTERMUTECTCALLS.out.versions)
    ch_versions = ch_versions.mix(GETPILEUPSUMMARIES.out.versions)
    ch_versions = ch_versions.mix(GATHERPILEUPSUMMARIES.out.versions)
    ch_versions = ch_versions.mix(LEARNREADORIENTATIONMODEL.out.versions)
    ch_versions = ch_versions.mix(MERGEMUTECTSTATS.out.versions)
    ch_versions = ch_versions.mix(MUTECT2.out.versions)

    emit:
    vcf   // channel: [ val(meta), [ vcf ] ]
    stats // channel: [ val(meta), [ stats ] ]

    artifact_priors     = LEARNREADORIENTATIONMODEL.out.artifactprior    // channel: [ val(meta), [ artifactprior ] ]

    pileup_table  // channel: [ val(meta), [ table ] ]

    contamination_table = CALCULATECONTAMINATION.out.contamination  // channel: [ val(meta), [ contamination ] ]
    segmentation_table  = CALCULATECONTAMINATION.out.segmentation   // channel: [ val(meta), [ segmentation ] ]

    vcf_filtered        = FILTERMUTECTCALLS.out.vcf.map{ meta, vcf -> [[
                                                                        id:             meta.sample,
                                                                        num_intervals:  meta.num_intervals,
                                                                        patient:        meta.patient,
                                                                        sample:         meta.sample,
                                                                        sex:            meta.sex,
                                                                        status:         meta.status,
                                                                        variantcaller:  "mutect2"
                                                                        ]
                                                                        , vcf] } // channel: [ val(meta), [ vcf ] ]
    filtered_index      = FILTERMUTECTCALLS.out.tbi                 // channel: [ val(meta), [ tbi ] ]
    filtered_stats      = FILTERMUTECTCALLS.out.stats               // channel: [ val(meta), [ stats ] ]

    versions            = ch_versions                               // channel: [ versions.yml ]
}
