//
// Manta germline variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { MANTA_GERMLINE } from '../../../modules/nf-core/manta/germline/main'

// Seems to be the consensus on upstream modules implementation too
workflow BAM_VARIANT_CALLING_GERMLINE_MANTA {
    take:
    cram          // channel: [mandatory] [ meta, cram, crai ]
    fasta         // channel: [mandatory] [ meta, fasta ]
    fasta_fai     // channel: [mandatory] [ meta, fasta_fai ]
    intervals     // channel: [mandatory] [ interval.bed.gz, interval.bed.gz.tbi] or [ [], []] if no intervals; intervals file contains all intervals

    main:
    versions = Channel.empty()

    // Combine cram and intervals, account for 0 intervals
    cram_intervals = cram.combine(intervals).map{ it ->
        bed_gz = it.size() > 3 ? it[3] : []
        bed_tbi = it.size() > 3 ? it[4] : []

        [it[0], it[1], it[2], bed_gz, bed_tbi]
    }

    MANTA_GERMLINE(cram_intervals, fasta, fasta_fai, [])

    small_indels_vcf = MANTA_GERMLINE.out.candidate_small_indels_vcf
    sv_vcf = MANTA_GERMLINE.out.candidate_sv_vcf
    diploid_sv_vcf = MANTA_GERMLINE.out.diploid_sv_vcf

    // Only diploid SV should get annotated
    // add variantcaller to meta map
    vcf = diploid_sv_vcf.map{ meta, vcf -> [ meta + [ variantcaller:'manta' ], vcf ] }

    versions = versions.mix(MANTA_GERMLINE.out.versions)

    emit:
    vcf

    versions
}
