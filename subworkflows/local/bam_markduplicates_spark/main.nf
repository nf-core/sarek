//
// MARKDUPLICATES SPARK AND QC after mapping
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { CRAM_QC_MOSDEPTH_SAMTOOLS                               } from '../cram_qc_mosdepth_samtools/main'
include { GATK4_ESTIMATELIBRARYCOMPLEXITY                         } from '../../../modules/nf-core/gatk4/estimatelibrarycomplexity/main'
include { GATK4SPARK_MARKDUPLICATES                               } from '../../../modules/nf-core/gatk4spark/markduplicates/main'
include { SAMTOOLS_INDEX                  as INDEX_MARKDUPLICATES } from '../../../modules/nf-core/samtools/index/main'

workflow BAM_MARKDUPLICATES_SPARK {
    take:
    bam                           // channel: [mandatory] meta, bam
    dict                          // channel: [mandatory] dict
    fasta                         // channel: [mandatory] fasta
    fasta_fai                     // channel: [mandatory] fasta_fai
    intervals_bed_combined        // channel: [optional]  intervals_bed

    main:
    reports = channel.empty()

    // RUN MARKUPDUPLICATES SPARK
    GATK4SPARK_MARKDUPLICATES(bam, fasta.map{ _meta, fasta_ -> [ fasta_ ] }, fasta_fai.map{ _meta, fasta_fai_ -> [ fasta_fai_ ] }, dict.map{ _meta, dict_ -> [ dict_ ] })

    // Index output (BAM or CRAM depending on ext.prefix)
    INDEX_MARKDUPLICATES(GATK4SPARK_MARKDUPLICATES.out.output)

    // Unified alignment output — join with the appropriate index
    alignment = GATK4SPARK_MARKDUPLICATES.out.output
        .join(INDEX_MARKDUPLICATES.out.index, failOnDuplicate: true, failOnMismatch: true)

    // QC on alignment
    CRAM_QC_MOSDEPTH_SAMTOOLS(alignment, fasta, fasta_fai, intervals_bed_combined)

    // When running Marduplicates spark, and saving reports
    GATK4_ESTIMATELIBRARYCOMPLEXITY(bam, fasta.map{ _meta, fasta_ -> [ fasta_ ] }, fasta_fai.map{ _meta, fasta_fai_ -> [ fasta_fai_ ] }, dict.map{ _meta, dict_ -> [ dict_ ] })

    // Gather all reports generated
    reports = reports.mix(GATK4_ESTIMATELIBRARYCOMPLEXITY.out.metrics)
    reports = reports.mix(CRAM_QC_MOSDEPTH_SAMTOOLS.out.reports)


    emit:
    alignment   // channel: [ meta, file, index ] — BAM or CRAM
    reports
}
