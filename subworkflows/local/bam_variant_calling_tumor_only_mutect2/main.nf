//
// Run GATK mutect2 in tumor only mode, getepileupsummaries, calculatecontamination and filtermutectcalls
//

include { GATK4_CALCULATECONTAMINATION    } from '../../../modules/nf-core/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS         } from '../../../modules/nf-core/gatk4/filtermutectcalls/main'
include { GATK4_GATHERPILEUPSUMMARIES     } from '../../../modules/nf-core/gatk4/gatherpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES        } from '../../../modules/nf-core/gatk4/getpileupsummaries/main'
include { GATK4_LEARNREADORIENTATIONMODEL } from '../../../modules/nf-core/gatk4/learnreadorientationmodel/main'
include { GATK4_MERGEMUTECTSTATS          } from '../../../modules/nf-core/gatk4/mergemutectstats/main'
include { GATK4_MERGEVCFS                 } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MUTECT2                   } from '../../../modules/nf-core/gatk4/mutect2/main'

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

    germline_resource_pileup     = germline_resource_tbi ? germline_resource : Channel.empty()
    germline_resource_pileup_tbi = germline_resource_tbi ?: Channel.empty()

    // Combine input and intervals for spread and gather strategy
    input_intervals = input.combine(intervals)
        // Move num_intervals to meta map and reorganize channel for GATK4_MUTECT2 module
        .map{ meta, input, index, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], input, index, intervals ] }

    if (joint_mutect2) {
        // Perform variant calling using mutect2 module in tumor single mode
        // Group cram files by patient
        patient_crams = input.map{ meta, t_cram, t_crai -> [ meta - meta.subMap('sample') + [id:meta.patient], t_cram, t_crai ] }.groupTuple()
        // Add intervals for scatter-gather scaling
        patient_cram_intervals = patient_crams.combine(intervals)
        // Move num_intervals to meta map and reorganize channel for GATK4_MUTECT2 module
            .map{ meta, t_cram, t_crai, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], t_cram, t_crai, intervals ] }
        GATK4_MUTECT2(patient_cram_intervals, fasta, fai, dict, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi)
    }
    else {
        // Perform variant calling using mutect2 module in tumor single mode
        GATK4_MUTECT2(input_intervals, fasta, fai, dict, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi)
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

    // Figuring out if there is one or more stats(s) from the same sample
    stats_branch = GATK4_MUTECT2.out.stats.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more f1r2(s) from the same sample
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

    // Mix intervals and no_intervals channels together
    vcf = Channel.empty().mix(GATK4_MERGEVCFS.out.vcf, vcf_branch.no_intervals)
    tbi = Channel.empty().mix(GATK4_MERGEVCFS.out.tbi, tbi_branch.no_intervals)
    stats = Channel.empty().mix(GATK4_MERGEMUTECTSTATS.out.stats, stats_branch.no_intervals)
    f1r2 = Channel.empty().mix(f1r2_to_merge, f1r2_branch.no_intervals)

    // Generate artifactpriors using learnreadorientationmodel on the f1r2 output of mutect2
    GATK4_LEARNREADORIENTATIONMODEL(f1r2)

    // Generate pileup summary table using getepileupsummaries
    GATK4_GETPILEUPSUMMARIES(input_intervals, fasta, fai, dict, germline_resource_pileup, germline_resource_pileup_tbi)

    // Figuring out if there is one or more table(s) from the same sample
    pileup_table_branch = GATK4_GETPILEUPSUMMARIES.out.table.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    pileup_table_to_merge = pileup_table_branch.intervals.map{ meta, table -> [ groupKey(meta, meta.num_intervals), table ] }.groupTuple()

    GATK4_GATHERPILEUPSUMMARIES(pileup_table_to_merge, dict.map{ meta, dict -> [ dict ] })

    // Mix intervals and no_intervals channels together
    pileup_table = Channel.empty().mix(GATK4_GATHERPILEUPSUMMARIES.out.table, pileup_table_branch.no_intervals)

    // Contamination and segmentation tables created using calculatecontamination on the pileup summary table
    GATK4_CALCULATECONTAMINATION(pileup_table.map{ meta, table -> [ meta, table, [] ] })

    // Reduce the meta to only patient name if joint_mutect2 otherwise keep regular ID
    contamination_table = joint_mutect2 ? GATK4_CALCULATECONTAMINATION.out.contamination.map{ meta, cont -> [ meta - meta.subMap('sample') + [id: meta.patient], cont]}.groupTuple() : GATK4_CALCULATECONTAMINATION.out.contamination
    segmentation_table  = joint_mutect2 ? GATK4_CALCULATECONTAMINATION.out.segmentation.map{  meta, seg  -> [ meta - meta.subMap('sample') + [id: meta.patient], seg]}.groupTuple()  : GATK4_CALCULATECONTAMINATION.out.segmentation

    // Mutect2 calls filtered by filtermutectcalls using the contamination and segmentation tables
    vcf_to_filter = vcf.join(tbi, failOnDuplicate: true, failOnMismatch: true)
        .join(stats, failOnDuplicate: true, failOnMismatch: true)
        .join(GATK4_LEARNREADORIENTATIONMODEL.out.artifactprior, failOnDuplicate: true, failOnMismatch: true)
        .join(segmentation_table, failOnDuplicate: true, failOnMismatch: true)
        .join(contamination_table, failOnDuplicate: true, failOnMismatch: true)
        .map{ meta, vcf, tbi, stats, artifactprior, seg, cont -> [ meta, vcf, tbi, stats, artifactprior, seg, cont, [] ] }

    GATK4_FILTERMUTECTCALLS(vcf_to_filter, fasta, fai, dict)

    vcf_filtered = GATK4_FILTERMUTECTCALLS.out.vcf
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'mutect2' ], vcf ] }

    versions = versions.mix(GATK4_CALCULATECONTAMINATION.out.versions)
    versions = versions.mix(GATK4_FILTERMUTECTCALLS.out.versions)
    versions = versions.mix(GATK4_GATHERPILEUPSUMMARIES.out.versions)
    versions = versions.mix(GATK4_GETPILEUPSUMMARIES.out.versions)
    versions = versions.mix(GATK4_LEARNREADORIENTATIONMODEL.out.versions)
    versions = versions.mix(GATK4_MERGEMUTECTSTATS.out.versions)
    versions = versions.mix(GATK4_MERGEVCFS.out.versions)
    versions = versions.mix(GATK4_MUTECT2.out.versions)

    emit:
    vcf   // channel: [ meta, vcf ]
    tbi   // channel: [ meta, tbi ]
    stats // channel: [ meta, stats ]

    vcf_filtered                                       // channel: [ meta, vcf ]
    index_filtered = GATK4_FILTERMUTECTCALLS.out.tbi   // channel: [ meta, tbi ]
    stats_filtered = GATK4_FILTERMUTECTCALLS.out.stats // channel: [ meta, stats ]

    artifact_priors = GATK4_LEARNREADORIENTATIONMODEL.out.artifactprior // channel: [ meta, artifactprior ]

    contamination_table // channel: [ meta, contamination ]
    pileup_table        // channel: [ meta, table ]
    segmentation_table  // channel: [ meta, segmentation ]

    versions // channel: [ versions.yml ]
}
