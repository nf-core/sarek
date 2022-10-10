//
// MARKDUPLICATES AND/OR QC after mapping
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BAM_TO_CRAM                            } from '../../bam_to_cram'
include { GATK4_ESTIMATELIBRARYCOMPLEXITY        } from '../../../../modules/nf-core/gatk4/estimatelibrarycomplexity/main'
include { GATK4_MARKDUPLICATES_SPARK             } from '../../../../modules/nf-core/gatk4/markduplicatesspark/main'
include { SAMTOOLS_INDEX as INDEX_MARKDUPLICATES } from '../../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_CONVERT as SAMTOOLS_CRAMTOBAM } from '../../../../modules/nf-core/samtools/convert/main'

workflow MARKDUPLICATES_SPARK {
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

    SAMTOOLS_CRAMTOBAM(cram_markduplicates, fasta, fasta_fai)

    // Convert Markupduplicates spark bam output to cram when running bamqc and/or deeptools
    BAM_TO_CRAM(Channel.empty(), cram_markduplicates, fasta, fasta_fai, intervals_bed_combined)

    // When running Marduplicates spark, and saving reports
    GATK4_ESTIMATELIBRARYCOMPLEXITY(bam, fasta, fasta_fai, dict)

    // Other reports done either with BAM_TO_CRAM subworkflow
    // or CRAM_QC subworkflow

    // Gather all reports generated
    qc_reports = qc_reports.mix(GATK4_ESTIMATELIBRARYCOMPLEXITY.out.metrics,
                                BAM_TO_CRAM.out.qc)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(GATK4_ESTIMATELIBRARYCOMPLEXITY.out.versions.first())
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES_SPARK.out.versions.first())
    ch_versions = ch_versions.mix(INDEX_MARKDUPLICATES.out.versions.first())
    ch_versions = ch_versions.mix(BAM_TO_CRAM.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_CRAMTOBAM.out.versions)

    emit:
        cram     = cram_markduplicates
        qc       = qc_reports

        versions = ch_versions // channel: [ versions.yml ]
}
