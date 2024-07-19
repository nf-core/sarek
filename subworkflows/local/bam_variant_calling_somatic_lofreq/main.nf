include { LOFREQ_SOMATIC  as  LOFREQSOMATIC } from '../../../modules/nf-core/lofreq/somatic/main.nf'
include { GATK4_MERGEVCFS as MERGE_LOFREQ_INDELS } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_LOFREQ_SNVS   } from '../../../modules/nf-core/gatk4/mergevcfs/main'

workflow BAM_VARIANT_CALLING_SOMATIC_LOFREQ {
    take:
    input     // channel: [mandatory] [ meta, normal_cram, normal_crai, tumor_cram, tumor_crai ]
    dict          // channel: [optional]  [ meta, dict ]
    fasta     // channel: [mandatory] [ fasta ]
    fai       // channel: [mandatory] [ fasta_fai ]
    intervals // channel: [mandatory] [ intervals ]

    main:
    versions = Channel.empty()

    input_intervals = input.combine(intervals)
        // Sort channel elements for LOFREQ_SOMATIC module
        .map {meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals, num_intervals -> [meta + [num_intervals:num_intervals], tumor_cram, tumor_crai, normal_cram, normal_crai, intervals]}
    
    LOFREQSOMATIC(input_intervals, fasta, fai)

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_indels = LOFREQSOMATIC.out.vcf_indels.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_snvs = LOFREQSOMATIC.out.vcf_snvs.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    vcf_indels_to_merge = vcf_indels.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()
    vcf_snvs_to_merge = vcf_snvs.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()

    MERGE_LOFREQ_INDELS(vcf_indels_to_merge, dict)
    MERGE_LOFREQ_SNVS(vcf_snvs_to_merge, dict)

    // Mix intervals and no_intervals channels together
    vcf = Channel.empty().mix(MERGE_LOFREQ_INDELS.out.vcf, MERGE_LOFREQ_SNVS.out.vcf, vcf_indels.no_intervals, vcf_snvs.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'lofreq' ], vcf ] }
    
    versions = versions.mix(LOFREQSOMATIC.out.versions)

    emit:
    vcf
    versions
}
