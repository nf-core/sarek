//
// dysgu variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { DYSGU } from '../../../modules/nf-core/dysgu/main'

// Seems to be the consensus on upstream modules implementation too
workflow BAM_VARIANT_CALLING_GERMLINE_DYSGU {
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

    DYSGU(cram_intervals, fasta, fasta_fai, [])

    
    dysgu_vcf = DYSGU.out.vcf

    // Only dysgu SV should get annotated
    // add variantcaller to meta map
    vcf = dysgu_vcf.map{ meta, vcf -> [ meta + [ variantcaller:'dysgu' ], vcf ] }

    versions = versions.mix(DYSGU.out.versions)

    emit:
    vcf

    versions
}
