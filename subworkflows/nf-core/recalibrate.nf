/*
========================================================================================
    RECALIBRATE
========================================================================================
*/

params.applybqsr_options      = [:]
params.applybqsr_spark_options = [:]
params.merge_cram_options      = [:]
params.qualimap_bamqc_options = [:]
params.samtools_index_options = [:]
params.samtools_stats_options = [:]

include { GATK4_APPLYBQSR as APPLYBQSR }             from '../../modules/nf-core/software/gatk4/applybqsr/main' addParams(options: params.applybqsr_options)
include { GATK4_APPLYBQSR_SPARK as APPLYBQSR_SPARK } from '../../modules/nf-core/software/gatk4/applybqsrspark/main' addParams(options: params.applybqsr_spark_options)
include { QUALIMAP_BAMQC }                           from '../../modules/nf-core/software/qualimap/bamqc/main'  addParams(options: params.qualimap_bamqc_options)
include { SAMTOOLS_INDEX }                           from '../../modules/nf-core/software/samtools/index/main'  addParams(options: params.samtools_index_options)
include { SAMTOOLS_MERGE_CRAM }                      from '../../modules/nf-core/software/samtools/merge_cram/main'  addParams(options: params.merge_cram_options)
include { SAMTOOLS_STATS }                           from '../../modules/nf-core/software/samtools/stats/main'  addParams(options: params.samtools_stats_options)

workflow RECALIBRATE {
    take:
        use_gatk_spark      //   value: [mandatory] use gatk spark
        skip_bamqc     // boolean: true/false
        skip_samtools  // boolean: true/false
        cram           // channel: [mandatory] cram
        dict           // channel: [mandatory] dict
        fai            // channel: [mandatory] fai
        fasta          // channel: [mandatory] fasta
        intervals      // channel: [mandatory] intervals
        target_bed     // channel: [optional]  target_bed

    main:

    cram_recalibrated_index = Channel.empty()
    cram_recalibrated       = Channel.empty()
    cram_reports            = Channel.empty()

    cram.combine(intervals).map{ meta, cram, crai, recal, intervals ->
        new_meta = meta.clone()
        new_meta.id = meta.sample + "_" + intervals.baseName
        [new_meta, cram, crai, recal, intervals]
    }.set{cram_intervals}

    if(use_gatk_spark){
        APPLYBQSR_SPARK(cram_intervals, fasta, fai, dict)
        cram_applybqsr = APPLYBQSR_SPARK.out.cram
    }else{
        APPLYBQSR(cram_intervals, fasta, fai, dict)
        cram_applybqsr = APPLYBQSR.out.cram
    }

    // STEP 4.5: MERGING AND INDEXING THE RECALIBRATED BAM FILES
    if (params.no_intervals) {
        cram_recalibrated = cram_applybqsr
    } else {
        cram_applybqsr.map{ meta, cram ->
            meta.id = meta.sample
            [meta, cram]
        }.groupTuple().set{cram_recalibrated_interval}

        SAMTOOLS_MERGE_CRAM(cram_recalibrated_interval)
        cram_recalibrated = SAMTOOLS_MERGE_CRAM.out.merged_cram

        SAMTOOLS_INDEX(cram_recalibrated)
        cram_recalibrated_index = cram_recalibrated.join(SAMTOOLS_INDEX.out.crai)

        qualimap_bamqc = Channel.empty()
        samtools_stats = Channel.empty()

        if (!skip_bamqc) {
            //TODO: Ugh BamQC again not liking crams, until we have a tool that can handle crams, not much choice but piping with samtools
            // QUALIMAP_BAMQC(bam_recalibrated, target_bed, params.target_bed)
            // qualimap_bamqc = QUALIMAP_BAMQC.out
        }

        if (!skip_samtools) {
            SAMTOOLS_STATS(cram_recalibrated_index, fasta)
            samtools_stats = SAMTOOLS_STATS.out.stats
        }
        cram_reports = samtools_stats.mix(qualimap_bamqc)
    }

    emit:
        cram = cram_recalibrated_index
        qc  = cram_reports
}
