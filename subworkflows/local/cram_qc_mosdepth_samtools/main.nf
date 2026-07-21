//
// QC on CRAM
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAMTOOLS_STATS } from '../../../modules/nf-core/samtools/stats/main'
include { MOSDEPTH       } from '../../../modules/nf-core/mosdepth/main'

workflow CRAM_QC_MOSDEPTH_SAMTOOLS {
    take:
    cram      // channel: [mandatory] [ meta, cram, crai ]
    fasta     // channel: [mandatory] [ meta, fasta ]
    fasta_fai // channel: [mandatory] [ meta, fai ]
    intervals

    main:
    versions = channel.empty()
    reports = channel.empty()

    // Reports run on cram
    SAMTOOLS_STATS(cram, fasta.combine(fasta_fai).map { meta, fasta_, _fai_meta, fai -> [ meta, fasta_, fai ] }.collect())

    MOSDEPTH(cram.combine(intervals.map { _meta, bed -> [bed ?: []] }), fasta, [])

    // Gather all reports generated
    reports = reports.mix(SAMTOOLS_STATS.out.stats)
    reports = reports.mix(MOSDEPTH.out.global_txt)
    reports = reports.mix(MOSDEPTH.out.regions_txt)

    // Gather versions of all tools used

    emit:
    reports
    versions // channel: [ versions.yml ]
}
