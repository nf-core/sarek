include { LOFREQ_CALLPARALLEL as LOFREQ       } from '../../../modules/nf-core/lofreq/callparallel/main.nf'
include { GATK4_MERGEVCFS     as MERGE_LOFREQ } from '../../../modules/nf-core/gatk4/mergevcfs/main.nf'

workflow BAM_VARIANT_CALLING_TUMOR_ONLY_LOFREQ {
    take:
    input     // channel: [mandatory] [ meta, tumor_cram, tumor_crai ]
    fasta     // channel: [mandatory] [ fasta ]
    fai       // channel: [mandatory] [ fasta_fai ]
    intervals // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ]
    dict      // channel: /path/to/reference/fasta/dictionary

    main:
    versions = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    input_intervals = input.combine(intervals)
        // Move num_intervals to meta map
        .map {meta, tumor_cram, tumor_crai, intervals, num_intervals -> [meta + [ num_intervals:num_intervals ], tumor_cram, tumor_crai, intervals]}

    LOFREQ(input_intervals, fasta, fai) // Call variants with LoFreq

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_branch = LOFREQ.out.vcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more tbi(s) from the same sample
    tbi_branch = LOFREQ.out.tbi.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    vcf_to_merge = vcf_branch.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ] }.groupTuple()

    MERGE_LOFREQ(vcf_to_merge, dict)

    // Mix intervals and no_intervals channels together
    // Remove unnecessary metadata
    vcf   = Channel.empty().mix(MERGE_LOFREQ.out.vcf, vcf_branch.no_intervals).map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'lofreq' ], vcf ] }

    versions = versions.mix(MERGE_LOFREQ.out.versions)
    versions = versions.mix(LOFREQ.out.versions)

    emit:
    vcf
    versions
}
