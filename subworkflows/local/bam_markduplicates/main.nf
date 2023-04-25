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
    bam                    // channel: [mandatory] [ meta, bam ]
    fasta                  // channel: [mandatory] [ fasta ]
    fasta_fai              // channel: [mandatory] [ fasta_fai ]
    intervals_bed_combined // channel: [optional]  [ intervals_bed ]

    main:
    versions = Channel.empty()
    reports  = Channel.empty()

    // RUN MARKUPDUPLICATES
    GATK4_MARKDUPLICATES(bam, fasta, fasta_fai)

    // Index cram
    INDEX_MARKDUPLICATES(GATK4_MARKDUPLICATES.out.cram)

    // Join with the crai file
    cram = GATK4_MARKDUPLICATES.out.cram.join(INDEX_MARKDUPLICATES.out.crai, failOnDuplicate: true, failOnMismatch: true)

    // QC on CRAM
    CRAM_QC_MOSDEPTH_SAMTOOLS(cram, fasta, intervals_bed_combined)

    // Gather all reports generated
    reports = reports.mix(GATK4_MARKDUPLICATES.out.metrics)
    reports = reports.mix(CRAM_QC_MOSDEPTH_SAMTOOLS.out.reports)

    // Gather versions of all tools used
    versions = versions.mix(GATK4_MARKDUPLICATES.out.versions)
    versions = versions.mix(INDEX_MARKDUPLICATES.out.versions)
    versions = versions.mix(CRAM_QC_MOSDEPTH_SAMTOOLS.out.versions)

    emit:
    cram
    reports

    versions    // channel: [ versions.yml ]
}
