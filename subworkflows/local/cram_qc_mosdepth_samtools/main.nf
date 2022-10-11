//
// QC on CRAM
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAMTOOLS_STATS     } from '../../../modules/nf-core/samtools/stats/main'
include { MOSDEPTH           } from '../../../modules/nf-core/mosdepth/main'

workflow CRAM_QC_MOSDEPTH_SAMTOOLS {
    take:
        cram                          // channel: [mandatory] meta, cram, crai
        fasta                         // channel: [mandatory] fasta
        fasta_fai                     // channel: [mandatory] fasta_fai
        intervals_bed_combined

    main:
    ch_versions = Channel.empty()
    qc_reports  = Channel.empty()

    // Reports run on cram
    SAMTOOLS_STATS(cram, fasta)
    MOSDEPTH(cram, intervals_bed_combined, fasta)

    // Gather all reports generated
    qc_reports = qc_reports.mix(SAMTOOLS_STATS.out.stats)
    qc_reports = qc_reports.mix(MOSDEPTH.out.global_txt,
                                MOSDEPTH.out.regions_txt)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())

    emit:
        qc       = qc_reports

        versions = ch_versions // channel: [ versions.yml ]
}
