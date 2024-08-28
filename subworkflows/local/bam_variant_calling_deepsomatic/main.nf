//
// DEEPSOMATIC variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { DEEPSOMATIC                               } from '../../../modules/nf-core/deepsomatic/main'
include { GATK4_MERGEVCFS as MERGE_DEEPSOMATIC_GVCF } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_DEEPSOMATIC_VCF  } from '../../../modules/nf-core/gatk4/mergevcfs/main'

// Deepvariant: https://github.com/google/deepvariant/issues/510
workflow BAM_VARIANT_CALLING_DEEPSOMATIC {
    take:
    cram_normal_tumor     // channel: [mandatory] [ meta, cram_normal, crai_normal, cram_tumor, crai_tumor ]
    dict          // channel: [optional]  [ meta, dict ]
    fasta         // channel: [mandatory] [ meta, fasta ]
    fasta_fai     // channel: [mandatory] [ meta, fasta_fai ]
    intervals     // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals

    main:
    versions = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    // Move num_intervals to meta map
    cram_normal_tumor_intervals = cram_normal_tumor.combine(intervals)
        // Move num_intervals to meta map
        .map{ meta, cram_n, crai_n, cram_t, crai_t, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], cram_n, crai_n, cram_t, crai_t ]}

    // Convert [intervals, num_intervals] to [meta, intervals] with an empty meta
    intervals_only = intervals.map { intervals, num_intervals -> [[], intervals]}

    DEEPSOMATIC(cram_normal_tumor_intervals, intervals_only, fasta, fasta_fai, [ [ id:'null' ], [] ])

    // // Figuring out if there is one or more vcf(s) from the same sample
    vcf_out = DEEPSOMATIC.out.vcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more gvcf(s) from the same sample
    gvcf_out = DEEPSOMATIC.out.gvcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // // Only when using intervals
    gvcf_to_merge = gvcf_out.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()
    vcf_to_merge = vcf_out.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()
    
    MERGE_DEEPSOMATIC_GVCF(gvcf_to_merge, dict)
    MERGE_DEEPSOMATIC_VCF(vcf_to_merge, dict)
    
    gvcf = Channel.empty()
    // Mix intervals and no_intervals channels together
    gvcf = Channel.empty().mix(MERGE_DEEPSOMATIC_GVCF.out.vcf, gvcf_out.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'deepsomatic' ], vcf ] }

    vcf = Channel.empty()
    // Mix intervals and no_intervals channels together
    vcf = Channel.empty().mix(MERGE_DEEPSOMATIC_VCF.out.vcf, vcf_out.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'deepsomatic' ], vcf ] }

    versions = versions.mix(DEEPSOMATIC.out.versions)
    // versions = versions.mix(MERGE_DEEPSOMATIC_GVCF.out.versions)
    versions = versions.mix(MERGE_DEEPSOMATIC_VCF.out.versions)

    emit:
    gvcf
    vcf
    versions
}
