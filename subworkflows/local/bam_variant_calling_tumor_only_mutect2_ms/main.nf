//
// Run GATK mutect2 in multi sample mode for tumor only samples, getpileupsummaries, calculatecontamination and filtermutectcalls
//

include { GATK4_MUTECT2                   as MUTECT2                   } from '../../../modules/nf-core/gatk4/mutect2/main'
include { GATK4_MERGEVCFS                 as MERGE_MUTECT2             } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEMUTECTSTATS          as MERGEMUTECTSTATS          } from '../../../modules/nf-core/gatk4/mergemutectstats/main'
include { GATK4_LEARNREADORIENTATIONMODEL as LEARNREADORIENTATIONMODEL } from '../../../modules/nf-core/gatk4/learnreadorientationmodel/main'
include { GATK4_GETPILEUPSUMMARIES        as GETPILEUPSUMMARIES        } from '../../../modules/nf-core/gatk4/getpileupsummaries/main'
include { GATK4_GATHERPILEUPSUMMARIES     as GATHERPILEUPSUMMARIES     } from '../../../modules/nf-core/gatk4/gatherpileupsummaries/main'
include { GATK4_CALCULATECONTAMINATION    as CALCULATECONTAMINATION    } from '../../../modules/nf-core/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS         as FILTERMUTECTCALLS         } from '../../../modules/nf-core/gatk4/filtermutectcalls/main'

workflow BAM_VARIANT_CALLING_TUMOR_ONLY_MUTECT2_MS {
    take:
    input                     // channel: [ val(meta), [ input ], [ input_index ] ]
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

    // Combine input and intervals for spread and gather strategy
    input_intervals = input.combine(intervals)
        // Move num_intervals to meta map and reorganize channel for MUTECT2 module
        .map{ meta, input, index, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], input, index, intervals ] }
    // Remove intervals, group tumor cram, crai files by patient.
    ch_tumor_cram = input_intervals.map{meta, t_cram, t_crai, intervals -> [[id: meta.patient, patient: meta.patient, sex: meta.sex, num_intervals: meta.num_intervals], t_cram, t_crai]}.unique().groupTuple()
    // Add intervals back
    ch_pt_intervals = input_intervals.map{meta, t_cram, t_crai, intervals -> [[id: meta.patient, patient: meta.patient, sex: meta.sex, num_intervals: meta.num_intervals], intervals]}.unique().groupTuple()
    ch_to_intervals = ch_tumor_cram.join(ch_pt_intervals).transpose(by : [3])
    //
    //Perform variant calling using mutect2 module in tumor single mode.
    //
    MUTECT2(ch_to_intervals,
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
        }.set{ mutect2_vcf_branch }

    MUTECT2.out.tbi.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{ mutect2_tbi_branch }

    MUTECT2.out.stats.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{ mutect2_stats_branch }

    MUTECT2.out.f1r2.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{ mutect2_f1r2_branch }

    //Only when using intervals
    //Merge Mutect2 VCF

    MERGE_MUTECT2(
        mutect2_vcf_branch.intervals
        .map{ meta, vcf ->
            new_meta = [id: meta.patient, patient: meta.patient, sex: meta.sex, num_intervals: meta.num_intervals]

            [groupKey(new_meta, meta.num_intervals), vcf]
        }.groupTuple(),
        dict)

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
            new_meta = [id: meta.patient, patient: meta.patient, sex: meta.sex, num_intervals: meta.num_intervals]
            [groupKey(new_meta, meta.num_intervals), stats]
        }.groupTuple())

    mutect2_stats = Channel.empty().mix(
        MERGEMUTECTSTATS.out.stats,
        mutect2_stats_branch.no_intervals)

    //
    //Generate artifactpriors using learnreadorientationmodel on the f1r2 output of mutect2.
    //
    LEARNREADORIENTATIONMODEL(
        Channel.empty().mix(
            mutect2_f1r2_branch.intervals
            .map{ meta, f1r2 ->
                new_meta = [id: meta.patient, patient: meta.patient, sex: meta.sex, num_intervals: meta.num_intervals]
                [groupKey(new_meta, meta.num_intervals), f1r2]
            }.groupTuple(),
            mutect2_f1r2_branch.no_intervals))

    //
    //Generate pileup summary table using getepileupsummaries.
    //
    
    GETPILEUPSUMMARIES ( input_intervals , fasta, fai, dict, germline_resource_pileup , germline_resource_pileup_tbi )

    GETPILEUPSUMMARIES.out.table.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }set{ pileup_table_branch }

    //Merge Pileup Summaries
    GATHERPILEUPSUMMARIES(
        pileup_table_branch.intervals
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
        dict.map{ meta, dict -> [ dict ]})

    pileup_table = Channel.empty().mix(
        GATHERPILEUPSUMMARIES.out.table,
        pileup_table_branch.no_intervals)

    //
    // Contamination and segmentation tables created using calculatecontamination on the pileup summary table.
    //
    table_contamination = pileup_table.map{meta, table -> [meta, table, []]}
    CALCULATECONTAMINATION ( table_contamination )

    //
    // Mutect2 calls filtered by filtermutectcalls using the contamination and segmentation tables.
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
    
    ch_filtermutect_in = ch_filtermutect.map{ meta, vcf, tbi, stats, artifactprior, seg, cont -> [meta, vcf, tbi, stats, artifactprior, seg, cont, []] }

    FILTERMUTECTCALLS ( ch_filtermutect_in, fasta, fai, dict )

    mutect2_vcf_filtered = FILTERMUTECTCALLS.out.vcf.map{ meta, vcf -> [[id:meta.id, variantcaller:"mutect2_ms"], vcf]}
    mutect2_vcf_filtered_tbi = FILTERMUTECTCALLS.out.tbi.map{ meta, tbi -> [[id:meta.id, variantcaller:"mutect2_ms"], tbi]}
    mutect2_vcf_filtered_stats = FILTERMUTECTCALLS.out.stats.map{ meta, stats -> [[id:meta.id, variantcaller:"mutect2_ms"], stats]}

    ch_versions = ch_versions.mix(MERGE_MUTECT2.out.versions)
    ch_versions = ch_versions.mix(CALCULATECONTAMINATION.out.versions)
    ch_versions = ch_versions.mix(FILTERMUTECTCALLS.out.versions)
    ch_versions = ch_versions.mix(GETPILEUPSUMMARIES.out.versions)
    ch_versions = ch_versions.mix(GATHERPILEUPSUMMARIES.out.versions)
    ch_versions = ch_versions.mix(LEARNREADORIENTATIONMODEL.out.versions)
    ch_versions = ch_versions.mix(MERGEMUTECTSTATS.out.versions)
    ch_versions = ch_versions.mix(MUTECT2.out.versions)

    emit:
    mutect2_vcf         = mutect2_vcf                                  // channel: [ val(meta), [ vcf ] ]
    mutect2_stats       = mutect2_stats                                // channel: [ val(meta), [ stats ] ]

    artifact_priors     = LEARNREADORIENTATIONMODEL.out.artifactprior  // channel: [ val(meta), [ artifactprior ] ]

    pileup_table        = pileup_table                                 // channel: [ val(meta), [ table ] ]

    contamination_table = CALCULATECONTAMINATION.out.contamination     // channel: [ val(meta), [ contamination ] ]
    segmentation_table  = CALCULATECONTAMINATION.out.segmentation      // channel: [ val(meta), [ segmentation ] ]

    filtered_vcf        = mutect2_vcf_filtered                         // channel: [ val(meta), [ vcf ] ]
    filtered_index      = mutect2_vcf_filtered_tbi                     // channel: [ val(meta), [ tbi ] ]
    filtered_stats      = mutect2_vcf_filtered_stats                   // channel: [ val(meta), [ stats ] ]

    versions            = ch_versions                                  // channel: [ versions.yml ]
}
