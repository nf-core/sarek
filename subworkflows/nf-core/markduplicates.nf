//
// MARKDUPLICATES AND/OR QC after mapping
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { DEEPTOOLS_BAMCOVERAGE                                    } from '../../modules/local/deeptools/bamcoverage'
include { GATK4_ESTIMATELIBRARYCOMPLEXITY                          } from '../../modules/nf-core/modules/gatk4/estimatelibrarycomplexity/main'
include { GATK4_MARKDUPLICATES                                     } from '../../modules/nf-core/modules/gatk4/markduplicates/main'
include { GATK4_MARKDUPLICATES_SPARK                               } from '../../modules/local/gatk4/markduplicatesspark/main'
include { QUALIMAP_BAMQC                                           } from '../../modules/local/qualimap/bamqc/main'
include { SAMTOOLS_INDEX as INDEX_MARKDUPLICATES                   } from '../../modules/local/samtools/index/main'
include { SAMTOOLS_STATS                                           } from '../../modules/nf-core/modules/samtools/stats/main'
include { SAMTOOLS_VIEWINDEX as SAMTOOLS_BAM_TO_CRAM_DUPLICATES    } from '../../modules/local/samtools/viewindex/main'
include { SAMTOOLS_VIEWINDEX as SAMTOOLS_BAM_TO_CRAM_NO_DUPLICATES } from '../../modules/local/samtools/viewindex/main'
include { SAMTOOLS_VIEWINDEX as SAMTOOLS_BAM_TO_CRAM_SPARK         } from '../../modules/local/samtools/viewindex/main'

workflow MARKDUPLICATES {
    take:
        bam_mapped                    // channel: [mandatory, when running Markduplicates, else optional] meta, bam
        bam_indexed                   // channel: [mandatory, when skipping Markduplicates, else optional] meta, bam, bai
        dict                          // channel: [mandatory] dict
        fasta                         // channel: [mandatory] fasta
        fasta_fai                     // channel: [mandatory] fasta_fai
        intervals_combined_bed_gz_tbi // channel: [optional]  intervals_bed.gz, intervals_bed.gz.tbi

    main:
    ch_versions = Channel.empty()
    qc_reports  = Channel.empty()

    // When skipping Markduplicates converting bam input to cram
    SAMTOOLS_BAM_TO_CRAM_NO_DUPLICATES(bam_indexed, fasta, fasta_fai)

    // Run Markupduplicates spark
    // When running bamqc and/or deeptools output is bam, else cram
    GATK4_MARKDUPLICATES_SPARK(bam_mapped, fasta, fasta_fai, dict)
    INDEX_MARKDUPLICATES(GATK4_MARKDUPLICATES_SPARK.out.output)

    // Convert Markupduplicates spark bam output to cram when running bamqc and/or deeptools
    SAMTOOLS_BAM_TO_CRAM_SPARK(GATK4_MARKDUPLICATES_SPARK.out.output.join(INDEX_MARKDUPLICATES.out.bam_bai), fasta, fasta_fai)

    // Run Markupduplicates
    GATK4_MARKDUPLICATES(bam_mapped)
    // Convert output to cram
    SAMTOOLS_BAM_TO_CRAM_DUPLICATES(GATK4_MARKDUPLICATES.out.bam.join(GATK4_MARKDUPLICATES.out.bai), fasta, fasta_fai)

    // Only one of these channel is not empty:
    // - skipping Markduplicates
    // - running Markupduplicates spark with bam output
    // - running Markupduplicates spark with cram output
    // - running Markupduplicates
    cram_markduplicates = Channel.empty().mix(
        SAMTOOLS_BAM_TO_CRAM_NO_DUPLICATES.out.cram_crai,
        SAMTOOLS_BAM_TO_CRAM_SPARK.out.cram_crai,
        GATK4_MARKDUPLICATES_SPARK.out.output.join(INDEX_MARKDUPLICATES.out.cram_crai),
        SAMTOOLS_BAM_TO_CRAM_DUPLICATES.out.cram_crai)

    // When running Marduplicates spark, and saving reports
    GATK4_ESTIMATELIBRARYCOMPLEXITY(bam_mapped, fasta, fasta_fai, dict)
    // Reports on Marduplicates spark bam output or on bam input
    QUALIMAP_BAMQC(bam_indexed.mix(GATK4_MARKDUPLICATES.out.bam.join(GATK4_MARKDUPLICATES.out.bai)), intervals_combined_bed_gz_tbi)
    DEEPTOOLS_BAMCOVERAGE(bam_indexed.mix(GATK4_MARKDUPLICATES.out.bam.join(GATK4_MARKDUPLICATES.out.bai)))
    // Other reports run on cram
    SAMTOOLS_STATS(cram_markduplicates, fasta)

    // Gather all reports generated
    qc_reports = qc_reports.mix(DEEPTOOLS_BAMCOVERAGE.out.bigwig)
    qc_reports = qc_reports.mix(GATK4_ESTIMATELIBRARYCOMPLEXITY.out.metrics)
    qc_reports = qc_reports.mix(GATK4_MARKDUPLICATES.out.metrics)
    qc_reports = qc_reports.mix(QUALIMAP_BAMQC.out.results)
    qc_reports = qc_reports.mix(SAMTOOLS_STATS.out.stats)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions.first())
    ch_versions = ch_versions.mix(GATK4_ESTIMATELIBRARYCOMPLEXITY.out.versions.first())
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions.first())
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES_SPARK.out.versions.first())
    ch_versions = ch_versions.mix(INDEX_MARKDUPLICATES.out.versions.first())
    ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_BAM_TO_CRAM_DUPLICATES.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_BAM_TO_CRAM_NO_DUPLICATES.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_BAM_TO_CRAM_SPARK.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())

    emit:
        cram     = cram_markduplicates
        qc       = qc_reports

        versions = ch_versions // channel: [ versions.yml ]
}
