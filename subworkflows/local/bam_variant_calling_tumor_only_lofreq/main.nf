include { LOFREQ_CALLPARALLEL } from '../../../modules/nf-core/lofreq/callparallel/main.nf'
include { BCFTOOLS_SORT as SORT } from '../../../modules/nf-core/bcftools/sort/main'
include { TABIX_TABIX as TABIX } from '../../../modules/nf-core/tabix/tabix/main'


workflow BAM_VARIANT_CALLING_TUMOR_ONLY_LOFREQ {
    take:
    input     // channel: [mandatory] [ meta, tumor_cram, tumor_crai ]
    fasta     // channel: [mandatory] [ fasta ]
    fai       // channel: [mandatory] [ fasta_fai ]
    intervals // channel: [mandatory] [ intervals ]

    main:
    versions = Channel.empty()

    input_intervals = input.combine(intervals)
        // Sort channel elements for LOFREQ_SOMATIC module
        .map {meta, tumor_cram, tumor_crai, intervals, num_intervals -> [meta + [ num_intervals:num_intervals ], tumor_cram, tumor_crai, intervals]}

    LOFREQ_CALLPARALLEL(input_intervals, fasta, fai)

    TABIX(LOFREQ_CALLPARALLEL.out.vcf)

    lofreq = LOFREQ_CALLPARALLEL.out.vcf.combine(TABIX.out.tbi).map {meta, vcf, meta2, tbi -> [meta, [vcf, tbi]]}

    SORT(lofreq)

    vcf = Channel.empty().mix(SORT.out.vcf)
        .map{ meta, vcf -> [ meta + [ variantcaller:'lofreq' ], vcf ] }
    versions = versions.mix(LOFREQ_CALLPARALLEL.out.versions)

    emit:
    vcf
    versions
}
