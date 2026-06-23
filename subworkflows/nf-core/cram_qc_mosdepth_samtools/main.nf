//
// QC on CRAM
//

include { MOSDEPTH       } from '../../../modules/nf-core/mosdepth'
include { SAMTOOLS_STATS } from '../../../modules/nf-core/samtools/stats'

workflow CRAM_QC_MOSDEPTH_SAMTOOLS {
    take:
    cram      // channel: [mandatory] [ meta, cram, crai ]
    fasta     // channel: [mandatory] [ fasta ]
    intervals

    main:
    reports = Channel.empty()

    // Reports run on cram
    SAMTOOLS_STATS(cram, fasta)

    MOSDEPTH(cram.combine(intervals.map { meta, bed -> [bed ?: []] }), fasta)

    // Gather all reports generated
    reports = reports.mix(SAMTOOLS_STATS.out.stats)
    reports = reports.mix(MOSDEPTH.out.global_txt)
    reports = reports.mix(MOSDEPTH.out.regions_txt)

    emit:
    reports // channel: [ meta, report_file ]
}
