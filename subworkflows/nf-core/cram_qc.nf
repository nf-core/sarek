//
// QC on CRAM
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAMTOOLS_STATS     } from '../../modules/nf-core/modules/samtools/stats/main'
include { QUALIMAP_BAMQCCRAM } from '../../modules/nf-core/modules/qualimap/bamqccram/main'

workflow CRAM_QC {
    take:
        cram                          // channel: [mandatory] meta, cram, crai
        fasta                         // channel: [mandatory] fasta
        fasta_fai                     // channel: [mandatory] fasta_fai
        intervals_combined_bed_gz_tbi

    main:
    ch_versions = Channel.empty()
    qc_reports  = Channel.empty()

    // Reports run on cram
    SAMTOOLS_STATS(cram, fasta)
    QUALIMAP_BAMQCCRAM(INDEX_RECALIBRATE.out.cram_crai, intervals_combined_bed_gz_tbi, fasta, fasta_fai)

    // Gather all reports generated
    qc_reports = qc_reports.mix(SAMTOOLS_STATS.out.stats)
    qc_reports = qc_reports.mix(QUALIMAP_BAMQCCRAM.out.results)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(QUALIMAP_BAMQCCRAM.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())

    emit:
        qc       = qc_reports

        versions = ch_versions // channel: [ versions.yml ]
}
