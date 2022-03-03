//
// MARKDUPLICATES
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_MARKDUPLICATES } from '../../modules/nf-core/modules/gatk4/markduplicates/main'
include { BAM_TO_CRAM          } from './bam_to_cram'

workflow MARKDUPLICATES {
    take:
        bam                           // channel: [mandatory] meta, bam
        fasta                         // channel: [mandatory] fasta
        fasta_fai                     // channel: [mandatory] fasta_fai
        intervals_combined_bed_gz_tbi // channel: [optional]  intervals_bed.gz, intervals_bed.gz.tbi

    main:
    ch_versions = Channel.empty()
    qc_reports  = Channel.empty()

    // Run Markupduplicates
    GATK4_MARKDUPLICATES(bam)

    // Convert output to cram
    BAM_TO_CRAM(GATK4_MARKDUPLICATES.out.bam.join(GATK4_MARKDUPLICATES.out.bai), fasta, fasta_fai, intervals_combined_bed_gz_tbi)

    // Gather all reports generated
    qc_reports = qc_reports.mix(GATK4_MARKDUPLICATES.out.metrics)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions.first())
    ch_versions = ch_versions.mix(BAM_TO_CRAM.out.versions)

    emit:
        cram     = BAM_TO_CRAM.out.cram
        qc       = qc_reports

        versions = ch_versions // channel: [ versions.yml ]
}
