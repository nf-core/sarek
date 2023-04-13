include { MANTA_SOMATIC                                    } from '../../../modules/nf-core/manta/somatic/main'

workflow BAM_VARIANT_CALLING_SOMATIC_MANTA {
    take:
    cram          // channel: [mandatory] [ meta, cram1, crai1, cram2, crai2 ]
    dict          // channel: [optional]  [ meta, dict ]
    fasta         // channel: [mandatory] [ fasta ]
    fasta_fai     // channel: [mandatory] [ fasta_fai ]
    intervals     // channel: [mandatory] [ interval.bed.gz, interval.bed.gz.tbi, num_intervals ] or [ [], [], 0 ] if no intervals

    main:
    versions = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    cram_intervals = cram.combine(intervals)

    MANTA_SOMATIC(cram_intervals, fasta, fasta_fai)

    // Figuring out if there is one or more vcf(s) from the same sample
    candidate_small_indels_vcf = MANTA_SOMATIC.out.candidate_small_indels_vcf

    // Figuring out if there is one or more vcf(s) from the same sample
    candidate_small_indels_vcf_tbi = MANTA_SOMATIC.out.candidate_small_indels_vcf_tbi

    // Figuring out if there is one or more vcf(s) from the same sample
    candidate_sv_vcf = MANTA_SOMATIC.out.candidate_sv_vcf

    // Figuring out if there is one or more vcf(s) from the same sample
    diploid_sv_vcf = MANTA_SOMATIC.out.diploid_sv_vcf

    // Figuring out if there is one or more vcf(s) from the same sample
    somatic_sv_vcf = MANTA_SOMATIC.out.somatic_sv_vcf


    // Mix intervals and no_intervals channels together
    // Only diploid and somatic SV should get annotated
    vcf = Channel.empty().mix(diploid_sv_vcf, somatic_sv_vcf)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta + [ variantcaller:'manta' ], vcf ] }

    versions = versions.mix(MANTA_SOMATIC.out.versions)

    emit:
    candidate_small_indels_vcf
    candidate_small_indels_vcf_tbi
    vcf

    versions
}
