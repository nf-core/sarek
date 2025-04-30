//
// MANTA somatic variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { MANTA_SOMATIC } from '../../../modules/nf-core/manta/somatic/main'

workflow BAM_VARIANT_CALLING_SOMATIC_MANTA {
    take:
    cram          // channel: [mandatory] [ meta, cram1, crai1, cram2, crai2 ]
    fasta         // channel: [mandatory] [ meta, fasta ]
    fasta_fai     // channel: [mandatory] [ meta, fasta_fai ]
    intervals     // channel: [mandatory] [ interval.bed.gz, interval.bed.gz.tbi ] or [ [], [] ] if no intervals

    main:
    versions = Channel.empty()

    // Combine cram and intervals, account for 0 intervals
    cram_intervals = cram.combine(intervals).map{ it ->
        bed_gz = it.size() > 5 ? it[5] : []
        bed_tbi = it.size() > 5 ? it[6] : []

        [it[0], it[1], it[2], it[3], it[4], bed_gz, bed_tbi]
    }

    MANTA_SOMATIC(cram_intervals, fasta, fasta_fai, [])

    candidate_small_indels_vcf = MANTA_SOMATIC.out.candidate_small_indels_vcf
    candidate_small_indels_vcf_tbi = MANTA_SOMATIC.out.candidate_small_indels_vcf_tbi
    candidate_sv_vcf = MANTA_SOMATIC.out.candidate_sv_vcf
    diploid_sv_vcf = MANTA_SOMATIC.out.diploid_sv_vcf
    somatic_sv_vcf = MANTA_SOMATIC.out.somatic_sv_vcf

    // Only diploid and somatic SV should get annotated
    // add variantcaller to meta map
    vcf = Channel.empty().mix(diploid_sv_vcf, somatic_sv_vcf).map{ meta, vcf -> [ meta + [ variantcaller:'manta' ], vcf ] }

    versions = versions.mix(MANTA_SOMATIC.out.versions)

    emit:
    candidate_small_indels_vcf
    candidate_small_indels_vcf_tbi
    vcf

    versions
}
