include { MANTA_SOMATIC                                    } from '../../../modules/nf-core/manta/somatic/main'

workflow BAM_VARIANT_CALLING_SOMATIC_MANTA {
    take:
    cram          // channel: [mandatory] [ meta, cram1, crai1, cram2, crai2 ]
    dict          // channel: [optional]  [ meta, dict ]
    fasta         // channel: [mandatory] [ fasta ]
    fasta_fai     // channel: [mandatory] [ fasta_fai ]
    intervals     // channel: [mandatory] [ interval.bed.gz, interval.bed.gz.tbi ] or [ [], [] ] if no intervals

    main:
    versions = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    cram_intervals = cram.combine(intervals)

    MANTA_SOMATIC(cram_intervals, fasta, fasta_fai)

    candidate_small_indels_vcf = MANTA_SOMATIC.out.candidate_small_indels_vcf
    candidate_small_indels_vcf_tbi = MANTA_SOMATIC.out.candidate_small_indels_vcf_tbi
    candidate_sv_vcf = MANTA_SOMATIC.out.candidate_sv_vcf
    diploid_sv_vcf = MANTA_SOMATIC.out.diploid_sv_vcf
    somatic_sv_vcf = MANTA_SOMATIC.out.somatic_sv_vcf

    // Only diploid and somatic SV should get annotated
    vcf = Channel.empty().mix(diploid_sv_vcf, somatic_sv_vcf)
        // add variantcaller to meta map
        .map{ meta, vcf -> [ meta + [ variantcaller:'manta' ], vcf ] }

    versions = versions.mix(MANTA_SOMATIC.out.versions)

    emit:
    candidate_small_indels_vcf
    candidate_small_indels_vcf_tbi
    vcf

    versions
}
