include { MANTA_TUMORONLY                                  } from '../../../modules/nf-core/manta/tumoronly/main'

// Seems to be the consensus on upstream modules implementation too
workflow BAM_VARIANT_CALLING_TUMOR_ONLY_MANTA {
    take:
    cram          // channel: [mandatory] [ meta, cram, crai ]
    dict          // channel: [optional]  [ meta, dict ]
    fasta         // channel: [mandatory] [ fasta ]
    fasta_fai     // channel: [mandatory] [ fasta_fai ]
    intervals     // channel: [mandatory] [ interval.bed.gz, interval.bed.gz.tbi, num_intervals ] or [ [], [], 0 ] if no intervals

    main:
    versions = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    cram_intervals = cram.combine(intervals)

    MANTA_TUMORONLY(cram_intervals, fasta, fasta_fai)

    // Figuring out if there is one or more vcf(s) from the same sample
    small_indels_vcf = MANTA_TUMORONLY.out.candidate_small_indels_vcf

    // Figuring out if there is one or more vcf(s) from the same sample
    candidate_sv_vcf = MANTA_TUMORONLY.out.candidate_sv_vcf

    // Figuring out if there is one or more vcf(s) from the same sample
    tumor_sv_vcf = MANTA_TUMORONLY.out.tumor_sv_vcf

    // Mix intervals and no_intervals channels together
    // Only tumor sv should get annotated
    vcf = tumor_sv_vcf
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta + [ variantcaller:'manta' ], vcf ] }

    versions = versions.mix(MANTA_TUMORONLY.out.versions)

    emit:
    vcf

    versions
}
