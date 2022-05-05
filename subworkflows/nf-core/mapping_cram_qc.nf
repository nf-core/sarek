//
// CRAM and optionnal QC for mapped reads provided as CRAM
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { DEEPTOOLS_BAMCOVERAGE } from '../../modules/nf-core/modules/deeptools/bamcoverage/main'
include { QUALIMAP_BAMQCCRAM } from '../../modules/nf-core/modules/qualimap/bamqccram/main'

workflow MAPPING_CRAM_QC {
    take:
        cram_indexed                  // channel: [manadtory] meta, cram, crai
        fasta                         // channel: [mandatory] fasta
        fasta_fai                     // channel: [mandatory] fai
        intervals_combined_bed_gz_tbi // channel: [optional]  intervals_bed.gz, intervals_bed.gz.tbi

    main:
    ch_versions = Channel.empty()
    qc_reports  = Channel.empty()

    // Reports on cram input
    DEEPTOOLS_BAMCOVERAGE(cram_indexed)
    QUALIMAP_BAMQCCRAM(cram_indexed, intervals_combined_bed_gz_tbi, fasta, fasta_fai)

    // Other reports run on cram

    // Gather all reports generated
    qc_reports = qc_reports.mix(DEEPTOOLS_BAMCOVERAGE.out.bigwig)
    qc_reports = qc_reports.mix(QUALIMAP_BAMQCCRAM.out.results)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions.first())
    ch_versions = ch_versions.mix(QUALIMAP_BAMQCCRAM.out.versions.first())

    emit:
        qc       = qc_reports

        versions = ch_versions // channel: [ versions.yml ]
}
