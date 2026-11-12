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
    skip_tools // string: --skip_tools

    main:
    reports = channel.empty()

    // Reports run on cram
    // `--skip_tools samtools` to skip samtools stats
    if (!(skip_tools && skip_tools.split(',').contains('samtools'))) {
        SAMTOOLS_STATS(cram, fasta.combine(fasta_fai).map { meta, fasta_, _fai_meta, fai -> [ meta, fasta_, fai ] }.collect())
        reports = reports.mix(SAMTOOLS_STATS.out.stats)
    }

    // `--skip_tools mosdepth` to skip mosdepth
    if (!(skip_tools && skip_tools.split(',').contains('mosdepth'))) {
        MOSDEPTH(cram.combine(intervals.map { _meta, bed -> [bed ?: []] }), fasta, [])
        reports = reports.mix(MOSDEPTH.out.global_txt)
        reports = reports.mix(MOSDEPTH.out.regions_txt)
    }

    emit:
    reports
}
