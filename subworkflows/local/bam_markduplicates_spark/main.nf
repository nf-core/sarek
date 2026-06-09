//
// MARKDUPLICATES SPARK AND QC after mapping
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { CRAM_QC_MOSDEPTH_SAMTOOLS       } from '../cram_qc_mosdepth_samtools/main'
include { GATK4_ESTIMATELIBRARYCOMPLEXITY } from '../../../modules/nf-core/gatk4/estimatelibrarycomplexity/main'
include { GATK4SPARK_MARKDUPLICATES       } from '../../../modules/nf-core/gatk4spark/markduplicates/main'

workflow BAM_MARKDUPLICATES_SPARK {
    take:
    bam                           // channel: [mandatory] meta, bam
    dict                          // channel: [mandatory] dict
    fasta                         // channel: [mandatory] fasta
    fasta_fai                     // channel: [mandatory] fasta_fai
    intervals_bed_combined        // channel: [optional]  intervals_bed

    main:
    versions = Channel.empty()
    reports = Channel.empty()

    // RUN MARKUPDUPLICATES SPARK
    // GATK4 MarkDuplicatesSpark indexes its output via --create-output-bam-index (default true);
    // module emits both .bai (BAM) and .crai (CRAM) on the bam_index channel.
    GATK4SPARK_MARKDUPLICATES(bam, fasta.map{ meta, fasta_ -> [ fasta_ ] }, fasta_fai.map{ meta, fasta_fai_ -> [ fasta_fai_ ] }, dict.map{ meta, dict_ -> [ dict_ ] })

    // Unified alignment output — BAM or CRAM depending on ext.prefix
    alignment = GATK4SPARK_MARKDUPLICATES.out.output
        .join(GATK4SPARK_MARKDUPLICATES.out.bam_index, failOnDuplicate: true, failOnMismatch: true)

    // QC on alignment
    CRAM_QC_MOSDEPTH_SAMTOOLS(alignment, fasta, intervals_bed_combined)

    // When running Marduplicates spark, and saving reports
    GATK4_ESTIMATELIBRARYCOMPLEXITY(bam, fasta.map{ meta, fasta_ -> [ fasta_ ] }, fasta_fai.map{ meta, fasta_fai_ -> [ fasta_fai_ ] }, dict.map{ meta, dict_ -> [ dict_ ] })

    // Gather all reports generated
    reports = reports.mix(GATK4_ESTIMATELIBRARYCOMPLEXITY.out.metrics)
    reports = reports.mix(CRAM_QC_MOSDEPTH_SAMTOOLS.out.reports)

    // Gather versions of all tools used
    versions = versions.mix(GATK4_ESTIMATELIBRARYCOMPLEXITY.out.versions)
    versions = versions.mix(GATK4SPARK_MARKDUPLICATES.out.versions)
    versions = versions.mix(CRAM_QC_MOSDEPTH_SAMTOOLS.out.versions)

    emit:
    alignment   // channel: [ meta, file, index ] — BAM or CRAM
    reports

    versions // channel: [ versions.yml ]
}
