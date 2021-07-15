/*
========================================================================================
    RECALIBRATE
========================================================================================
*/

params.applybqsr_options      = [:]
params.merge_bam_options      = [:]
params.qualimap_bamqc_options = [:]
params.samtools_index_options = [:]
params.samtools_stats_options = [:]

include { GATK4_APPLYBQSR as APPLYBQSR } from '../../modules/nf-core/software/gatk4/applybqsr/main' addParams(options: params.applybqsr_options)
include { QUALIMAP_BAMQC }               from '../../modules/nf-core/software/qualimap/bamqc/main'  addParams(options: params.qualimap_bamqc_options)
include { SAMTOOLS_INDEX }               from '../../modules/nf-core/software/samtools/index/main'  addParams(options: params.samtools_index_options)
include { SAMTOOLS_MERGE }               from '../../modules/nf-core/software/samtools/merge/main'  addParams(options: params.merge_bam_options)
include { SAMTOOLS_STATS }               from '../../modules/nf-core/software/samtools/stats/main'  addParams(options: params.samtools_stats_options)

workflow RECALIBRATE {
    take:
        skip_bamqc     // boolean: true/false
        skip_samtools  // boolean: true/false
        bam            // channel: [mandatory] bam
        dict           // channel: [mandatory] dict
        fai            // channel: [mandatory] fai
        fasta          // channel: [mandatory] fasta
        intervals      // channel: [mandatory] intervals
        target_bed     // channel: [optional]  target_bed

    main:

    bam_recalibrated_index = Channel.empty()
    bam_recalibrated       = Channel.empty()
    bam_reports            = Channel.empty()

    bam.combine(intervals).map{ meta, bam, bai, recal, intervals ->
        new_meta = meta.clone()
        new_meta.id = meta.sample + "_" + intervals.baseName
        [new_meta, bam, bai, recal, intervals]
    }.set{bam_intervals}

    APPLYBQSR(bam_intervals, fasta, fai, dict)

    // STEP 4.5: MERGING AND INDEXING THE RECALIBRATED BAM FILES
    if (params.no_intervals) {
        bam_recalibrated = APPLYBQSR.out.bam
    } else {
        APPLYBQSR.out.bam.map{ meta, bam ->
            meta.id = meta.sample
            [meta, bam]
        }.groupTuple().set{bam_recalibrated_interval}

        SAMTOOLS_MERGE(bam_recalibrated_interval)
        bam_recalibrated = SAMTOOLS_MERGE.out.merged_bam

        SAMTOOLS_INDEX(bam_recalibrated)
        bam_recalibrated_index = bam_recalibrated.join(SAMTOOLS_INDEX.out.bai)

        qualimap_bamqc = Channel.empty()
        samtools_stats = Channel.empty()

        if (!skip_bamqc) {
            QUALIMAP_BAMQC(bam_recalibrated, target_bed, params.target_bed)
            qualimap_bamqc = QUALIMAP_BAMQC.out
        }

        if (!skip_samtools) {
            SAMTOOLS_STATS(bam_recalibrated_index)
            samtools_stats = SAMTOOLS_STATS.out.stats
        }
        bam_reports = samtools_stats.mix(qualimap_bamqc)
    }

    emit:
        bam = bam_recalibrated_index
        qc  = bam_reports
}
