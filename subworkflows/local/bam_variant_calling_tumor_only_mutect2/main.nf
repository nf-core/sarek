//
// GATK MUTECT2 in tumor only mode: getepileupsummaries, calculatecontamination and filtermutectcalls
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

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

    //If no germline resource is provided, then create an empty channel to avoid GetPileupsummaries from being run
    germline_resource_pileup     = germline_resource_tbi ? germline_resource : Channel.empty()
    germline_resource_pileup_tbi = germline_resource_tbi ?: Channel.empty()

    // Combine input and intervals for spread and gather strategy
    input_intervals = input.combine(intervals)
        // Move num_intervals to meta map and reorganize channel for MUTECT2 module
        .map{ meta, input, index, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], input, index, intervals ] }

    if (joint_mutect2) {
        // Perform variant calling using mutect2 module in tumor single mode
        // Group cram files by patient
        input_joint = input
            .map{ meta, input, index -> [ meta - meta.subMap('sample') + [ id:meta.patient ], input, index ] }
            .groupTuple()

        // Add intervals for scatter-gather scaling
        input_joint_intervals = input_joint.combine(intervals)
        // Move num_intervals to meta map and reorganize channel for MUTECT2 module
            .map{ meta, cram, crai, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], cram, crai, intervals ] }
        MUTECT2(input_joint_intervals, fasta, fai, dict, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi)
    }
    else {
        // Perform variant calling using mutect2 module in tumor single mode
        MUTECT2(input_intervals, fasta, fai, dict, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi)
    }

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_branch = MUTECT2.out.vcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more tbi(s) from the same sample
    tbi_branch = MUTECT2.out.tbi.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more stats(s) from the same sample
    stats_branch = MUTECT2.out.stats.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more f1r2(s) from the same sample
    f1r2_branch = MUTECT2.out.f1r2.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    vcf_to_merge = vcf_branch.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ] }.groupTuple()
    stats_to_merge = stats_branch.intervals.map{ meta, stats -> [ groupKey(meta, meta.num_intervals), stats ] }.groupTuple()
    f1r2_to_merge = f1r2_branch.intervals.map{ meta, f1r2 -> [ groupKey(meta, meta.num_intervals), f1r2 ] }.groupTuple()

    MERGE_MUTECT2(vcf_to_merge, dict)
    MERGEMUTECTSTATS(stats_to_merge)

    // Mix intervals and no_intervals channels together
    // Remove unnecessary metadata
    vcf   = Channel.empty().mix(MERGE_MUTECT2.out.vcf, vcf_branch.no_intervals).map{ meta, vcf -> [ meta - meta.subMap('num_intervals'), vcf ] }
    tbi   = Channel.empty().mix(MERGE_MUTECT2.out.tbi, tbi_branch.no_intervals).map{ meta, tbi -> [ meta - meta.subMap('num_intervals'), tbi ] }
    stats = Channel.empty().mix(MERGEMUTECTSTATS.out.stats, stats_branch.no_intervals).map{ meta, stats -> [ meta - meta.subMap('num_intervals'), stats ] }
    f1r2  = Channel.empty().mix(f1r2_to_merge, f1r2_branch.no_intervals).map{ meta, f1r2 -> [ meta - meta.subMap('num_intervals'), f1r2 ] }

    // Generate artifactpriors using learnreadorientationmodel on the f1r2 output of mutect2
    LEARNREADORIENTATIONMODEL(f1r2)

    pileup_input = input_intervals.map{ meta, cram, crai, intervals -> [ meta + [ id:meta.sample ], cram, crai, intervals] }.unique()

    // Generate pileup summary table using getepileupsummaries
    GETPILEUPSUMMARIES(pileup_input, fasta, fai, dict, germline_resource_pileup, germline_resource_pileup_tbi)

    // Figuring out if there is one or more table(s) from the same sample
    pileup_table_branch = GETPILEUPSUMMARIES.out.table.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    pileup_table_to_merge = pileup_table_branch.intervals.map{ meta, table -> [ groupKey(meta, meta.num_intervals), table ] }.groupTuple()

    GATHERPILEUPSUMMARIES(pileup_table_to_merge, dict.map{ meta, dict -> [ dict ] })

    // Mix intervals and no_intervals channels together
    pileup_table = Channel.empty().mix(GATHERPILEUPSUMMARIES.out.table, pileup_table_branch.no_intervals).map{meta, table -> [ meta - meta.subMap('num_intervals') + [id:meta.sample], table ] }

    // Contamination and segmentation tables created using calculatecontamination on the pileup summary table
    CALCULATECONTAMINATION(pileup_table.map{ meta, table -> [ meta, table, [] ] })

    // Initialize empty channel: Contamination calculation is run on pileup table, pileup is not run if germline resource is not provided
    calculatecontamination_out_seg = Channel.empty()
    calculatecontamination_out_cont = Channel.empty()

    if (joint_mutect2) {
        // Group tables by samples
        calculatecontamination_out_seg = CALCULATECONTAMINATION.out.segmentation.map{ meta, seg -> [ meta - meta.subMap('sample', 'num_intervals') + [id:meta.patient], seg ] }.groupTuple()
        calculatecontamination_out_cont = CALCULATECONTAMINATION.out.contamination.map{ meta, cont -> [ meta - meta.subMap('sample', 'num_intervals') + [id:meta.patient], cont ] }.groupTuple()
    } else {
        // Regular single sample mode
        calculatecontamination_out_seg = CALCULATECONTAMINATION.out.segmentation.map{ meta, seg -> [ meta - meta.subMap('num_intervals'), seg ] }
        calculatecontamination_out_cont = CALCULATECONTAMINATION.out.contamination.map{ meta, cont -> [ meta - meta.subMap('num_intervals'), cont ] }
    }

    // Mutect2 calls filtered by filtermutectcalls using the contamination and segmentation tables
    vcf_to_filter = vcf.join(tbi, failOnDuplicate: true, failOnMismatch: true)
        .join(stats, failOnDuplicate: true, failOnMismatch: true)
        .join(LEARNREADORIENTATIONMODEL.out.artifactprior, failOnDuplicate: true, failOnMismatch: true)
        .join(calculatecontamination_out_seg)
        .join(calculatecontamination_out_cont)
        .map{ meta, vcf, tbi, stats, artifactprior, seg, cont -> [ meta, vcf, tbi, stats, artifactprior, seg, cont, [] ] }

    FILTERMUTECTCALLS(vcf_to_filter, fasta, fai, dict)

    vcf_filtered = FILTERMUTECTCALLS.out.vcf
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta + [ variantcaller:'mutect2' ], vcf ] }

    versions = versions.mix(MERGE_MUTECT2.out.versions)
    versions = versions.mix(CALCULATECONTAMINATION.out.versions)
    versions = versions.mix(FILTERMUTECTCALLS.out.versions)
    versions = versions.mix(GETPILEUPSUMMARIES.out.versions)
    versions = versions.mix(GATHERPILEUPSUMMARIES.out.versions)
    versions = versions.mix(LEARNREADORIENTATIONMODEL.out.versions)
    versions = versions.mix(MERGEMUTECTSTATS.out.versions)
    versions = versions.mix(MUTECT2.out.versions)

    emit:
    vcf   // channel: [ meta, vcf ]
    stats // channel: [ meta, stats ]

    vcf_filtered                                 // channel: [ meta, vcf ]
    index_filtered = FILTERMUTECTCALLS.out.tbi   // channel: [ meta, tbi ]
    stats_filtered = FILTERMUTECTCALLS.out.stats // channel: [ meta, stats ]

    artifact_priors = LEARNREADORIENTATIONMODEL.out.artifactprior    // channel: [ meta, artifactprior ]

    pileup_table  // channel: [ meta, table ]

    contamination_table = calculatecontamination_out_cont  // channel: [ meta, contamination ]
    segmentation_table  = calculatecontamination_out_seg   // channel: [ meta, segmentation ]

    versions // channel: [ versions.yml ]
}
