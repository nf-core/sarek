include { LOFREQ_CALLPARALLEL } from '../../../modules/nf-core/lofreq/callparallel/main.nf'
include { BEDTOOLS_SORT as SORT } from '../../../modules/nf-core/bedtools/sort'
include { BEDTOOLS_MERGE as MERGE } from '../../../modules/nf-core/bedtools/merge'



workflow BAM_VARIANT_CALLING_TUMOR_ONLY_LOFREQ {
    take:
    input     // channel: [mandatory] [ meta, tumor_cram, tumor_crai ]
    fasta     // channel: [mandatory] [ fasta ]
    fai       // channel: [mandatory] [ fasta_fai ]
    intervals // channel: [path] [ intervals ]

    main:
    versions = Channel.empty()

    //Modify the .bed file to avoid errors with duplicate variants.
    fasta_genome = fai.map {meta, fai -> [fai]}
    intervals_sort = intervals.map { bed, interval -> [ [ id:bed.baseName + "_unsorted" ], bed ] }
    SORT(intervals_sort, fasta_genome)
    interval_merge = SORT.out.sorted.map {meta, sorted -> [[id:sorted.baseName + "_unmerge"], sorted]}
    MERGE(interval_merge)
    intervals_final = intervals.combine(MERGE.out.bed).map {path, inter, meta, bed -> [bed, inter]}

    // Combine cram and intervals for spread and gather strategy
    input_intervals = input.combine(intervals_final)
        // Move num_intervals to meta map
        .map {meta, tumor_cram, tumor_crai, intervals, num_intervals -> [meta + [ num_intervals:num_intervals ], tumor_cram, tumor_crai, intervals]}

    LOFREQ_CALLPARALLEL(input_intervals, fasta, fai)

    vcf = Channel.empty().mix(LOFREQ_CALLPARALLEL.out.vcf)
        .map{ meta, vcf -> [ meta + [ variantcaller:'lofreq' ], vcf ] }
    versions = versions.mix(LOFREQ_CALLPARALLEL.out.versions)

    emit:
    vcf
    versions
}
