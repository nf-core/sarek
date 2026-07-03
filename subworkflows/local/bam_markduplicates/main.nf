//
// MARKDUPLICATES AND QC after mapping
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { CRAM_QC_MOSDEPTH_SAMTOOLS } from '../cram_qc_mosdepth_samtools/main'
include { GATK4_MARKDUPLICATES      } from '../../../modules/nf-core/gatk4/markduplicates/main'

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
    // --CREATE_INDEX true is set via ext.args when --save_output_as_bam, so the
    // module emits .bai inline; CRAM mode emits .crai via samtools post-conversion.
    GATK4_MARKDUPLICATES(bam, fasta.map{ meta, fasta_ -> [ fasta_ ] }, fasta_fai.map{ meta, fasta_fai_ -> [ fasta_fai_ ] })

    // Unified alignment output — BAM or CRAM depending on save_output_as_bam
    alignment = GATK4_MARKDUPLICATES.out.bam
        .join(GATK4_MARKDUPLICATES.out.bai, failOnDuplicate: true, failOnMismatch: true)
        .mix(GATK4_MARKDUPLICATES.out.cram
            .join(GATK4_MARKDUPLICATES.out.crai, failOnDuplicate: true, failOnMismatch: true))

    // QC on alignment
    CRAM_QC_MOSDEPTH_SAMTOOLS(alignment, fasta, intervals_bed_combined)

    // Gather all reports generated
    reports = reports.mix(GATK4_MARKDUPLICATES.out.metrics)
    reports = reports.mix(CRAM_QC_MOSDEPTH_SAMTOOLS.out.reports)

    // Gather versions of all tools used
    versions = versions.mix(GATK4_MARKDUPLICATES.out.versions)
    versions = versions.mix(CRAM_QC_MOSDEPTH_SAMTOOLS.out.versions)

    emit:
    alignment   // channel: [ meta, file, index ] — BAM or CRAM
    reports

    versions    // channel: [ versions.yml ]
}
