//
// MARKDUPLICATES AND/OR QC after mapping
//

params.estimatelibrarycomplexity_options = [:]
params.markduplicates_options            = [:]
params.markduplicatesspark_options       = [:]
params.merge_bam_options                 = [:]
params.qualimap_bamqc_options            = [:]
params.samtools_index_options            = [:]
params.samtools_stats_options            = [:]
params.samtools_view_options             = [:]

include { GATK4_ESTIMATELIBRARYCOMPLEXITY }             from '../../modules/local/gatk4/estimatelibrarycomplexity/main' addParams(options: params.estimatelibrarycomplexity_options)
include { GATK4_MARKDUPLICATES }                        from '../../modules/local/gatk4/markduplicates/main'            addParams(options: params.markduplicates_options)
include { GATK4_MARKDUPLICATES_SPARK }                  from '../../modules/local/gatk4/markduplicatesspark/main'       addParams(options: params.markduplicatesspark_options)
include { QUALIMAP_BAMQC }                              from '../../modules/local/qualimap/bamqc/main'                  addParams(options: params.qualimap_bamqc_options)
include { SAMTOOLS_INDEX }                              from '../../modules/local/samtools/index/main'                  addParams(options: params.samtools_index_options)
include { SAMTOOLS_STATS }                              from '../../modules/local/samtools/stats/main'                  addParams(options: params.samtools_stats_options)
include { SAMTOOLS_VIEW as SAMTOOLS_BAM_TO_CRAM }       from '../../modules/local/samtools/view/main'                   addParams(options: params.samtools_view_options)
include { SAMTOOLS_VIEW as SAMTOOLS_BAM_TO_CRAM_SPARK } from '../../modules/local/samtools/view/main'                   addParams(options: params.samtools_view_options)

workflow MARKDUPLICATES {
    take:
        bam_mapped          // channel: [mandatory] meta, bam
        bam_indexed         // channel: [mandatory] meta, bam, bai
        use_gatk_spark      //   value: [mandatory] use gatk spark
        save_metrics        //   value: [mandatory] save metrics
        fasta               // channel: [mandatory] fasta
        fai                 // channel: [mandatory] fai
        dict                // channel: [mandatory] dict
        skip_markduplicates // boolean: true/false
        skip_bamqc          // boolean: true/false
        skip_samtools       // boolean: true/false
        target_bed          // channel: [optional]  target_bed

    main:

    report_markduplicates = Channel.empty()

    if (skip_markduplicates) {
        bam_markduplicates = bam_indexed
        SAMTOOLS_BAM_TO_CRAM(bam_markduplicates, fasta, fai)
        cram_markduplicates = SAMTOOLS_BAM_TO_CRAM.out.cram_crai
    } else {
        if (use_gatk_spark) {
            //If BAMQC should be run on MD output, then don't use MDSpark to convert to cram, but use bam output instead
            if (!skip_bamqc) {
                GATK4_MARKDUPLICATES_SPARK(bam_mapped, fasta, fai, dict, "bam")
                SAMTOOLS_INDEX(GATK4_MARKDUPLICATES_SPARK.out.output)
                bam_markduplicates  = GATK4_MARKDUPLICATES_SPARK.out.output.join(SAMTOOLS_INDEX.out.bai)

                SAMTOOLS_BAM_TO_CRAM_SPARK(bam_markduplicates, fasta, fai)
                cram_markduplicates = SAMTOOLS_BAM_TO_CRAM_SPARK.out.cram_crai
            } else {
                GATK4_MARKDUPLICATES_SPARK(bam_mapped, fasta, fai, dict, "cram")
                SAMTOOLS_INDEX(GATK4_MARKDUPLICATES_SPARK.out.output)
                cram_markduplicates = GATK4_MARKDUPLICATES_SPARK.out.output.join(SAMTOOLS_INDEX.out.crai)
            }

            if (save_metrics) {
                GATK4_ESTIMATELIBRARYCOMPLEXITY(bam_mapped, fasta, fai, dict)
                report_markduplicates = GATK4_ESTIMATELIBRARYCOMPLEXITY.out.metrics
            }

        } else {
            GATK4_MARKDUPLICATES(bam_mapped)
            report_markduplicates = GATK4_MARKDUPLICATES.out.metrics
            bam_markduplicates    = GATK4_MARKDUPLICATES.out.bam_bai

            SAMTOOLS_BAM_TO_CRAM(bam_markduplicates, fasta, fai)
            cram_markduplicates = SAMTOOLS_BAM_TO_CRAM.out.cram_crai
        }
    }

    //If skip_markduplicates then QC tools are run on mapped bams,
    //if !skip_markduplicates, then QC tools are run on duplicate marked crams
    //After bamqc finishes, convert to cram for further analysis
    samtools_stats = Channel.empty()
    if (!skip_samtools) {
        SAMTOOLS_STATS(cram_markduplicates, fasta)
        samtools_stats = SAMTOOLS_STATS.out.stats
    }

    qualimap_bamqc = Channel.empty()
    if (!skip_bamqc) {
        QUALIMAP_BAMQC(bam_markduplicates, target_bed, params.target_bed)
        qualimap_bamqc = QUALIMAP_BAMQC.out.results
    }

    qc_reports = samtools_stats.mix(qualimap_bamqc)
    qc_reports = report_markduplicates.mix(qc_reports)

    emit:
        cram     = cram_markduplicates
        qc       = qc_reports
}
