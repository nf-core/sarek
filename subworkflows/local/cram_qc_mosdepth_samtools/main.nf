//
// QC on CRAM
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAMTOOLS_STATS     } from '../../../modules/nf-core/samtools/stats/main'
include { MOSDEPTH           } from '../../../modules/nf-core/mosdepth/main'

workflow CRAM_QC_MOSDEPTH_SAMTOOLS {
    take:
    cram                          // channel: [mandatory] [ meta, cram, crai ]
    fasta                         // channel: [mandatory] [ fasta ]
    intervals

    main:
    versions = Channel.empty()
    reports = Channel.empty()

    // Reports run on cram
    SAMTOOLS_STATS(cram, fasta)

    MOSDEPTH(cram.combine(intervals.map{ meta, bed -> [ bed?:[] ] }), fasta)

    // Gather all reports generated
    reports = reports.mix(SAMTOOLS_STATS.out.stats)
    reports = reports.mix(MOSDEPTH.out.global_txt)
    reports = reports.mix(MOSDEPTH.out.regions_txt)

    // Gather versions of all tools used
    versions = versions.mix(MOSDEPTH.out.versions)
    versions = versions.mix(SAMTOOLS_STATS.out.versions.first())

    emit:
    reports

    versions // channel: [ versions.yml ]
}
