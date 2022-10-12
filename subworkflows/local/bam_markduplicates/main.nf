//
// MARKDUPLICATES
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { CRAM_QC_MOSDEPTH_SAMTOOLS              } from '../cram_qc_mosdepth_samtools/main'
include { GATK4_MARKDUPLICATES                   } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { SAMTOOLS_INDEX as INDEX_MARKDUPLICATES } from '../../../modules/nf-core/samtools/index/main'

workflow BAM_MARKDUPLICATES {
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
    INDEX_MARKDUPLICATES(GATK4_MARKDUPLICATES.out.cram)

    cram_markduplicates = GATK4_MARKDUPLICATES.out.cram
        .join(INDEX_MARKDUPLICATES.out.crai)

    // Convert output to cram
    CRAM_QC_MOSDEPTH_SAMTOOLS(cram_markduplicates, fasta, fasta_fai, intervals_bed_combined)

    // Gather all reports generated
    qc_reports = qc_reports.mix(GATK4_MARKDUPLICATES.out.metrics,
                                CRAM_QC_MOSDEPTH_SAMTOOLS.out.qc)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions)
    ch_versions = ch_versions.mix(INDEX_MARKDUPLICATES.out.versions)
    ch_versions = ch_versions.mix(CRAM_QC_MOSDEPTH_SAMTOOLS.out.versions)

    emit:
        cram     = cram_markduplicates
        qc       = qc_reports

        versions = ch_versions // channel: [ versions.yml ]
}
