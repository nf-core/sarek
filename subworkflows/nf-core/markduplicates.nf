//
// MARKDUPLICATES AND/OR QC after mapping
//

include { GATK4_ESTIMATELIBRARYCOMPLEXITY                  } from '../../modules/nf-core/modules/gatk4/estimatelibrarycomplexity/main'
include { GATK4_MARKDUPLICATES                             } from '../../modules/nf-core/modules/gatk4/markduplicates/main'
include { GATK4_MARKDUPLICATES_SPARK                       } from '../../modules/local/gatk4/markduplicatesspark/main'
include { QUALIMAP_BAMQC                                   } from '../../modules/local/qualimap/bamqc/main'
include { SAMTOOLS_INDEX as INDEX_MARKDUPLICATES           } from '../../modules/local/samtools/index/main'
include { SAMTOOLS_STATS                                   } from '../../modules/nf-core/modules/samtools/stats/main'
include { SAMTOOLS_VIEWINDEX as SAMTOOLS_BAM_TO_CRAM       } from '../../modules/local/samtools/viewindex/main'
include { SAMTOOLS_VIEWINDEX as SAMTOOLS_BAM_TO_CRAM_SPARK } from '../../modules/local/samtools/viewindex/main'
include { DEEPTOOLS_BAMCOVERAGE                            } from '../../modules/local/deeptools/bamcoverage'

workflow MARKDUPLICATES {
    take:
        bam_mapped          // channel: [mandatory, if --skip_markdiplicate is false, else optional] meta, bam
        bam_indexed         // channel: [mandatory, if --skip_markduplicates is set, else optional] meta, bam, bai
        use_gatk_spark      //   value: [mandatory] use gatk spark
        save_metrics        //   value: [mandatory] save metrics
        dict                // channel: [mandatory] dict
        fasta               // channel: [mandatory] fasta
        fasta_fai           // channel: [mandatory] fasta_fai
        skip_markduplicates // boolean: true/false
        skip_bamqc          // boolean: true/false
        skip_samtools       // boolean: true/false
        intervals           // channel: [optional]  target_bed

    main:

    ch_versions           = Channel.empty()
    report_markduplicates = Channel.empty()

    if (skip_markduplicates) {
        bam_bai_markduplicates = bam_indexed
        SAMTOOLS_BAM_TO_CRAM(bam_bai_markduplicates, fasta, fasta_fai)
        cram_markduplicates = SAMTOOLS_BAM_TO_CRAM.out.cram_crai

        ch_versions = ch_versions.mix(SAMTOOLS_BAM_TO_CRAM.out.versions.first())
    } else {
        if (use_gatk_spark) {
            //If BAMQC should be run on MD output, then don't use MDSpark to convert to cram, but use bam output instead
            if (!skip_bamqc) {
                GATK4_MARKDUPLICATES_SPARK(bam_mapped, fasta, fasta_fai, dict, "bam")
                INDEX_MARKDUPLICATES(GATK4_MARKDUPLICATES_SPARK.out.output)
                bam_markduplicates = GATK4_MARKDUPLICATES_SPARK.out.output.join(INDEX_MARKDUPLICATES.out.bam_bai)

                SAMTOOLS_BAM_TO_CRAM_SPARK(bam_markduplicates, fasta, fasta_fai)
                cram_markduplicates = SAMTOOLS_BAM_TO_CRAM_SPARK.out.cram_crai

                ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES_SPARK.out.versions.first())
                ch_versions = ch_versions.mix(INDEX_MARKDUPLICATES.out.versions.first())
                ch_versions = ch_versions.mix(SAMTOOLS_BAM_TO_CRAM_SPARK.out.versions.first())
            } else {
                GATK4_MARKDUPLICATES_SPARK(bam_mapped, fasta, fasta_fai, dict, "cram")
                INDEX_MARKDUPLICATES(GATK4_MARKDUPLICATES_SPARK.out.output)
                cram_markduplicates = GATK4_MARKDUPLICATES_SPARK.out.output.join(INDEX_MARKDUPLICATES.out.crai)

                ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES_SPARK.out.versions.first())
                ch_versions = ch_versions.mix(INDEX_MARKDUPLICATES.out.versions.first())
            }

            if (save_metrics) {
                GATK4_ESTIMATELIBRARYCOMPLEXITY(bam_mapped, fasta, fasta_fai, dict)
                report_markduplicates = GATK4_ESTIMATELIBRARYCOMPLEXITY.out.metrics

                ch_versions = ch_versions.mix(GATK4_ESTIMATELIBRARYCOMPLEXITY.out.versions.first())
            }

        } else {
            GATK4_MARKDUPLICATES(bam_mapped)
            report_markduplicates  = GATK4_MARKDUPLICATES.out.metrics
            bam_markduplicates     = GATK4_MARKDUPLICATES.out.bam
            bai_markduplicates     = GATK4_MARKDUPLICATES.out.bai
            bam_bai_markduplicates = bam_markduplicates.join(bai_markduplicates)

            SAMTOOLS_BAM_TO_CRAM(bam_bai_markduplicates, fasta, fasta_fai)
            cram_markduplicates = SAMTOOLS_BAM_TO_CRAM.out.cram_crai

            ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions.first())
            ch_versions = ch_versions.mix(SAMTOOLS_BAM_TO_CRAM.out.versions.first())
        }
    }

    //If skip_markduplicates then QC tools are run on mapped bams,
    //if !skip_markduplicates, then QC tools are run on duplicate marked crams
    //After bamqc finishes, convert to cram for further analysis
    samtools_stats = Channel.empty()
    if (!skip_samtools) {
        SAMTOOLS_STATS(cram_markduplicates, fasta)
        samtools_stats = SAMTOOLS_STATS.out.stats

        ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())
    }

    qualimap_bamqc = Channel.empty()
    if (!skip_bamqc) {
        //TODO: intervals also with WGS data? Probably need a parameter if WGS for deepvariant tool, that would allow to check here too
        //TODO: error when no_intervals is set
        QUALIMAP_BAMQC(bam_bai_markduplicates, intervals)
        qualimap_bamqc = QUALIMAP_BAMQC.out.results

        ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions.first())
    }

    qc_reports = samtools_stats.mix(qualimap_bamqc)
    qc_reports = report_markduplicates.mix(qc_reports)

    emit:
        cram     = cram_markduplicates
        qc       = qc_reports

        versions = ch_versions // channel: [ versions.yml ]
}
