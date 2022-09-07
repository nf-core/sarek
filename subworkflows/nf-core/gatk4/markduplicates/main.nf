//
// MARKDUPLICATES
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_MARKDUPLICATES } from '../../../../modules/nf-core/modules/gatk4/markduplicates/main'
include { CRAM_QC              } from '../../cram_qc'
include { SAMTOOLS_INDEX       } from '../../../../modules/nf-core/modules/samtools/index/main'

workflow MARKDUPLICATES {
    take:
        bam                           // channel: [mandatory] meta, bam
        fasta                         // channel: [mandatory] fasta
        fasta_fai                     // channel: [mandatory] fasta_fai
        intervals_bed_combined        // channel: [optional]  intervals_bed

    main:
    ch_versions = Channel.empty()
    qc_reports  = Channel.empty()

    // Run Markupduplicates
    GATK4_MARKDUPLICATES(bam, fasta, fasta_fai)
    SAMTOOLS_INDEX(GATK4_MARKDUPLICATES.out.cram)

    cram_markduplicates = GATK4_MARKDUPLICATES.out.cram
        .join(SAMTOOLS_INDEX.out.crai)

    // Convert output to cram
    CRAM_QC(cram_markduplicates, fasta, fasta_fai, intervals_bed_combined)

    // Gather all reports generated
    qc_reports = qc_reports.mix(GATK4_MARKDUPLICATES.out.metrics,
                                CRAM_QC.out.qc)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions.first())
    ch_versions = ch_versions.mix(CRAM_QC.out.versions)

    emit:
        cram     = cram_markduplicates
        qc       = qc_reports

        versions = ch_versions // channel: [ versions.yml ]
}
