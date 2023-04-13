include { MANTA_GERMLINE                              } from '../../../modules/nf-core/manta/germline/main'

// Seems to be the consensus on upstream modules implementation too
workflow BAM_VARIANT_CALLING_GERMLINE_MANTA {
    take:
    cram          // channel: [mandatory] [ meta, cram, crai ]
    dict          // channel: [optional]  [ meta, dict ]
    fasta         // channel: [mandatory] [ fasta ]
    fasta_fai     // channel: [mandatory] [ fasta_fai ]
    intervals     // channel: [mandatory] [ interval.bed.gz, interval.bed.gz.tbi] or [ [], []] if no intervals; intervals file contains all intervals

    main:
    versions = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    cram_intervals = cram.combine(intervals)

    MANTA_GERMLINE(cram_intervals, fasta, fasta_fai)

    small_indels_vcf = MANTA_GERMLINE.out.candidate_small_indels_vcf
    sv_vcf = MANTA_GERMLINE.out.candidate_sv_vcf
    diploid_sv_vcf = MANTA_GERMLINE.out.diploid_sv_vcf

    // Only diploid SV should get annotated
    vcf = diploid_sv_vcf.
        // add variantcaller to meta map
        map{ meta, vcf -> [ meta + [ variantcaller:'manta' ], vcf ] }

    versions = versions.mix(MANTA_GERMLINE.out.versions)

    emit:
    vcf

    versions
}
