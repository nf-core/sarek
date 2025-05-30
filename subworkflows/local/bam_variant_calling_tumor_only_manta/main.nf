//
// MANTA single sample variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { MANTA_TUMORONLY } from '../../../modules/nf-core/manta/tumoronly/main'

// Seems to be the consensus on upstream modules implementation too
workflow BAM_VARIANT_CALLING_TUMOR_ONLY_MANTA {
    take:
    cram          // channel: [mandatory] [ meta, cram, crai ]
    fasta         // channel: [mandatory] [ meta, fasta ]
    fasta_fai     // channel: [mandatory] [ meta, fasta_fai ]
    intervals     // channel: [mandatory] [ interval.bed.gz, interval.bed.gz.tbi ] or [ [], [] ] if no intervals

    main:
    versions = Channel.empty()

    // Combine cram and intervals, account for 0 intervals
    cram_intervals = cram.combine(intervals).map{ it ->
        def bed_gz = it.size() > 3 ? it[3] : []
        def bed_tbi = it.size() > 3 ? it[4] : []

        [it[0], it[1], it[2], bed_gz, bed_tbi]
    }

    MANTA_TUMORONLY(cram_intervals, fasta, fasta_fai, [])

    small_indels_vcf = MANTA_TUMORONLY.out.candidate_small_indels_vcf
    candidate_sv_vcf = MANTA_TUMORONLY.out.candidate_sv_vcf
    tumor_sv_vcf = MANTA_TUMORONLY.out.tumor_sv_vcf

    // Only tumor sv should get annotated
    // add variantcaller to meta map
    vcf = tumor_sv_vcf.map{ meta, vcf -> [ meta + [ variantcaller:'manta' ], vcf ] }

    versions = versions.mix(MANTA_TUMORONLY.out.versions)

    emit:
    vcf

    versions
}
