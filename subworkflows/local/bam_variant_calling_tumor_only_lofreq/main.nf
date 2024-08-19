include { LOFREQ_CALLPARALLEL as LOFREQ          } from '../../../modules/nf-core/lofreq/callparallel/main.nf'

workflow BAM_VARIANT_CALLING_TUMOR_ONLY_LOFREQ {
    take:
    input     // channel: [mandatory] [ meta, tumor_cram, tumor_crai ]
    fasta     // channel: [mandatory] [ fasta ]
    fai       // channel: [mandatory] [ fasta_fai ]
    intervals // channel: [path] [ intervals ]

    main:
    versions = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    input_intervals = input.combine(intervals)
        // Move num_intervals to meta map
        .map {meta, tumor_cram, tumor_crai, intervals, num_intervals -> [meta + [ num_intervals:num_intervals ], tumor_cram, tumor_crai, intervals]}

    LOFREQ(input_intervals, fasta, fai) // Call variants with LoFreq

    vcf = Channel.empty().mix(LOFREQ.out.vcf)
        .map{ meta, vcf -> [ meta + [ variantcaller:'lofreq' ], vcf ] }

    versions = versions.mix(LOFREQ.out.versions)

    emit:
    vcf
    versions
}
