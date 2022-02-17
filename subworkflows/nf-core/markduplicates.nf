//
// MARKDUPLICATES AND/OR QC after mapping
//

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
        bam_mapped                      // channel: [mandatory, if no --skip_tools markduplicates, else optional] meta, bam
        bam_indexed                     // channel: [mandatory, if --skip_tools markduplicates, else optional] meta, bam, bai
        dict                            // channel: [mandatory] dict
        fasta                           // channel: [mandatory] fasta
        fasta_fai                       // channel: [mandatory] fasta_fai
        intervals_combined_bed_gz_tbi   // channel: [optional]  intervals_bed.gz, intervals_bed.gz.tbi

    main:

    ch_versions = Channel.empty()
    qc_reports  = Channel.empty()

    SAMTOOLS_BAM_TO_CRAM_NO_DUPLICATES(bam_indexed, fasta, fasta_fai)
    cram_markduplicates = SAMTOOLS_BAM_TO_CRAM_NO_DUPLICATES.out.cram_crai

    ch_versions = ch_versions.mix(SAMTOOLS_BAM_TO_CRAM_NO_DUPLICATES.out.versions.first())

    GATK4_MARKDUPLICATES_SPARK(bam_mapped, fasta, fasta_fai, dict)
    INDEX_MARKDUPLICATES(GATK4_MARKDUPLICATES_SPARK.out.output)

    //If BAMQC should be run on MD output, then don't use MDSpark to convert to cram, but use bam output instead
    bam_markduplicates  = GATK4_MARKDUPLICATES_SPARK.out.output.join(INDEX_MARKDUPLICATES.out.bam_bai)
    cram_markduplicates = GATK4_MARKDUPLICATES_SPARK.out.output.join(INDEX_MARKDUPLICATES.out.cram_crai)

    SAMTOOLS_BAM_TO_CRAM_SPARK(bam_markduplicates, fasta, fasta_fai)
    cram_markduplicates = SAMTOOLS_BAM_TO_CRAM_SPARK.out.cram_crai

    cram_markduplicates = GATK4_MARKDUPLICATES_SPARK.out.output.join(INDEX_MARKDUPLICATES.out.cram_crai)

    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES_SPARK.out.versions.first())
    ch_versions = ch_versions.mix(INDEX_MARKDUPLICATES.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_BAM_TO_CRAM_SPARK.out.versions.first())

    GATK4_MARKDUPLICATES(bam_mapped)
    SAMTOOLS_BAM_TO_CRAM_DUPLICATES(GATK4_MARKDUPLICATES.out.bam.join(GATK4_MARKDUPLICATES.out.bai), fasta, fasta_fai)

    cram_markduplicates = SAMTOOLS_BAM_TO_CRAM_DUPLICATES.out.cram_crai

    qc_reports = qc_reports.mix(GATK4_MARKDUPLICATES.out.metrics)

    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_BAM_TO_CRAM_DUPLICATES.out.versions.first())

    // GATK4_ESTIMATELIBRARYCOMPLEXITY is only run with SPARK, otherwise reports is from GATK4_MARKDUPLICATES
    // If --skip_tools markduplicates then QC tools are run on mapped bams,
    // If no --skip_tools markduplicates, then QC tools are run on duplicate marked crams
    // After bamqc finishes, convert to cram for further analysis

    GATK4_ESTIMATELIBRARYCOMPLEXITY(bam_mapped, fasta, fasta_fai, dict)
    SAMTOOLS_STATS(cram_markduplicates, fasta)
    QUALIMAP_BAMQC(GATK4_MARKDUPLICATES.out.bam.join(GATK4_MARKDUPLICATES.out.bai), intervals_combined_bed_gz_tbi)
    DEEPTOOLS_BAMCOVERAGE(GATK4_MARKDUPLICATES.out.bam.join(GATK4_MARKDUPLICATES.out.bai))

    qc_reports = qc_reports.mix(GATK4_ESTIMATELIBRARYCOMPLEXITY.out.metrics)
    qc_reports = qc_reports.mix(SAMTOOLS_STATS.out.stats)
    qc_reports = qc_reports.mix(QUALIMAP_BAMQC.out.results)
    qc_reports = qc_reports.mix(DEEPTOOLS_BAMCOVERAGE.out.bigwig)

    ch_versions = ch_versions.mix(GATK4_ESTIMATELIBRARYCOMPLEXITY.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions)
    ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions.first())

    emit:
        cram     = cram_markduplicates
        qc       = qc_reports

        versions = ch_versions // channel: [ versions.yml ]
}
