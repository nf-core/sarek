include { LOFREQ_CALLPARALLEL as LOFREQ          } from '../../../modules/nf-core/lofreq/callparallel/main.nf'
include { BEDTOOLS_SORT       as SORT_INTERVALS  } from '../../../modules/nf-core/bedtools/sort/main.nf'
include { BEDTOOLS_MERGE      as MERGE_INTERVALS } from '../../../modules/nf-core/bedtools/merge/main.nf'

workflow BAM_VARIANT_CALLING_TUMOR_ONLY_LOFREQ {
    take:
    input     // channel: [mandatory] [ meta, tumor_cram, tumor_crai ]
    fasta     // channel: [mandatory] [ fasta ]
    fai       // channel: [mandatory] [ fasta_fai ]
    intervals // channel: [path] [ intervals ]

    main:
    versions = Channel.empty()

    // Just making sure the intervals.bed file is properly forme to avoid error in lofreq.
    fasta_genome = fai.map {meta, fai -> [fai]}
    intervals_sort = intervals.map { bed, interval -> [ [ id:bed.baseName + "_unsorted_unmerged" ], bed ] }

    SORT_INTERVALS(intervals_sort, fasta_genome) // Sort the intervals

    interval_merge = SORT_INTERVALS.out.sorted.map {meta, sorted -> [[id:sorted.baseName + "_sorted_unmerged"], sorted]}

    MERGE_INTERVALS(interval_merge) // Merge overlapping intervals (if any)

    intervals_final = intervals.combine(MERGE_INTERVALS.out.bed).map {path, inter, meta, bed -> [bed, inter]}

    // Combine cram and intervals for spread and gather strategy
    input_intervals = input.combine(intervals_final)
        // Move num_intervals to meta map
        .map {meta, tumor_cram, tumor_crai, intervals, num_intervals -> [meta + [ num_intervals:num_intervals ], tumor_cram, tumor_crai, intervals]}

    LOFREQ(input_intervals, fasta, fai) // Call variants with LoFreq

    vcf = Channel.empty().mix(LOFREQ.out.vcf)
        .map{ meta, vcf -> [ meta + [ variantcaller:'lofreq' ], vcf ] }

    versions = versions.mix(LOFREQ.out.versions)
    versions = versions.mix(SORT_INTERVALS.out.versions)
    versions = versions.mix(MERGE_INTERVALS.out.versions)

    emit:
    vcf
    versions
}
