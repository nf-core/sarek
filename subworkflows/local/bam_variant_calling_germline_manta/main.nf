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
    versions = channel.empty()

    // Combine cram and intervals, account for 0 intervals
    cram_intervals = cram.combine(intervals).map{ combined ->
        def bed_gz = combined.size() > 3 ? combined[3] : []
        def bed_tbi = combined.size() > 3 ? combined[4] : []

        [combined[0], combined[1], combined[2], bed_gz, bed_tbi]
    }

    MANTA_GERMLINE(cram_intervals, fasta, fasta_fai, [])

    // add variantcaller to meta map
    candidate_small_indels_vcf     = MANTA_GERMLINE.out.candidate_small_indels_vcf.map{ meta, vcf -> [ meta + [ variantcaller:'manta' ], vcf ] }
    candidate_small_indels_vcf_tbi = MANTA_GERMLINE.out.candidate_small_indels_vcf_tbi.map{ meta, tbi -> [ meta + [ variantcaller:'manta' ], tbi ] }
    candidate_sv_vcf               = MANTA_GERMLINE.out.candidate_sv_vcf.map{ meta, vcf -> [ meta + [ variantcaller:'manta' ], vcf ] }
    candidate_sv_vcf_tbi           = MANTA_GERMLINE.out.candidate_sv_vcf_tbi.map{ meta, tbi -> [ meta + [ variantcaller:'manta' ], tbi ] }
    diploid_sv_vcf                 = MANTA_GERMLINE.out.diploid_sv_vcf.map{ meta, vcf -> [ meta + [ variantcaller:'manta' ], vcf ] }
    diploid_sv_vcf_tbi             = MANTA_GERMLINE.out.diploid_sv_vcf_tbi.map{ meta, tbi -> [ meta + [ variantcaller:'manta' ], tbi ] }


    emit:
    candidate_small_indels_vcf
    candidate_small_indels_vcf_tbi
    candidate_sv_vcf
    candidate_sv_vcf_tbi
    diploid_sv_vcf
    diploid_sv_vcf_tbi

    versions
}
