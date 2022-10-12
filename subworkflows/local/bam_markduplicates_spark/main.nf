//
// MARKDUPLICATES AND/OR QC after mapping
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { CRAM_QC_MOSDEPTH_SAMTOOLS              } from '../cram_qc_mosdepth_samtools/main'
include { GATK4_ESTIMATELIBRARYCOMPLEXITY        } from '../../../modules/nf-core/gatk4/estimatelibrarycomplexity/main'
include { GATK4_MARKDUPLICATES_SPARK             } from '../../../modules/nf-core/gatk4/markduplicatesspark/main'
include { SAMTOOLS_INDEX as INDEX_MARKDUPLICATES } from '../../../modules/nf-core/samtools/index/main'

workflow BAM_MARKDUPLICATES_SPARK {
    take:
        bam                           // channel: [mandatory] meta, bam
        dict                          // channel: [mandatory] dict
        fasta                         // channel: [mandatory] fasta
        fasta_fai                     // channel: [mandatory] fasta_fai
        intervals_bed_combined        // channel: [optional]  intervals_bed

    main:
    ch_versions = Channel.empty()
    qc_reports  = Channel.empty()

    // Run Markupduplicates spark
    // When running bamqc and/or deeptools output is bam, else cram
    GATK4_MARKDUPLICATES_SPARK(bam, fasta, fasta_fai, dict)
    INDEX_MARKDUPLICATES(GATK4_MARKDUPLICATES_SPARK.out.output)

    cram_markduplicates = GATK4_MARKDUPLICATES_SPARK.out.output
        .join(INDEX_MARKDUPLICATES.out.crai)

    // Convert Markupduplicates spark bam output to cram when running bamqc and/or deeptools
    CRAM_QC_MOSDEPTH_SAMTOOLS(cram_markduplicates, fasta, fasta_fai, intervals_bed_combined)

    // When running Marduplicates spark, and saving reports
    GATK4_ESTIMATELIBRARYCOMPLEXITY(bam, fasta, fasta_fai, dict)

    // Gather all reports generated
    qc_reports = qc_reports.mix(GATK4_ESTIMATELIBRARYCOMPLEXITY.out.metrics,
                                CRAM_QC_MOSDEPTH_SAMTOOLS.out.qc)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(GATK4_ESTIMATELIBRARYCOMPLEXITY.out.versions.first())
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES_SPARK.out.versions)
    ch_versions = ch_versions.mix(INDEX_MARKDUPLICATES.out.versions)
    ch_versions = ch_versions.mix(CRAM_QC_MOSDEPTH_SAMTOOLS.out.versions)

    emit:
        cram     = cram_markduplicates
        qc       = qc_reports

        versions = ch_versions // channel: [ versions.yml ]
}
