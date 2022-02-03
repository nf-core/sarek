//
// RECALIBRATE
//

include { GATK4_APPLYBQSR as APPLYBQSR              } from '../../modules/nf-core/modules/gatk4/applybqsr/main'
include { GATK4_APPLYBQSR_SPARK as APPLYBQSR_SPARK  } from '../../modules/local/gatk4/applybqsrspark/main'
include { QUALIMAP_BAMQC_CRAM                       } from '../../modules/local/qualimap/bamqccram/main'
include { SAMTOOLS_INDEX as INDEX_RECALIBRATE       } from '../../modules/local/samtools/index/main'
include { SAMTOOLS_MERGE_CRAM                       } from '../../modules/local/samtools/mergecram/main'
include { SAMTOOLS_STATS                            } from '../../modules/nf-core/modules/samtools/stats/main'

workflow RECALIBRATE {
    take:
        use_gatk_spark //   value: [mandatory] use gatk spark
        skip_bamqc     // boolean: true/false
        skip_samtools  // boolean: true/false
        cram           // channel: [mandatory] cram
        dict           // channel: [mandatory] dict
        fasta          // channel: [mandatory] fasta
        fasta_fai      // channel: [mandatory] fasta_fai
        intervals      // channel: [mandatory] intervals
        num_intervals
        no_intervals
        intervals_combined_bed_gz_tbi
        //target_bed     // channel: [optional]  target_bed

    main:

    ch_versions           = Channel.empty()
    cram_recalibrated_index = Channel.empty()
    cram_recalibrated       = Channel.empty()
    cram_reports            = Channel.empty()

    cram.combine(intervals).map{ meta, cram, crai, recal, intervals ->
        new_meta = meta.clone()
        new_meta.id = intervals.baseName != "no_intervals" ? meta.sample + "_" + intervals.baseName : meta.sample
        [new_meta, cram, crai, recal, intervals]
    }.set{cram_intervals}

    if(use_gatk_spark){
        APPLYBQSR_SPARK(cram_intervals, fasta, fasta_fai, dict)
        cram_applybqsr = APPLYBQSR_SPARK.out.cram
        ch_versions = ch_versions.mix(APPLYBQSR_SPARK.out.versions)
    }else{
        APPLYBQSR(cram_intervals, fasta, fasta_fai, dict)
        cram_applybqsr = APPLYBQSR.out.cram
        ch_versions = ch_versions.mix(APPLYBQSR.out.versions)
    }

    // STEP 4.5: MERGING AND INDEXING THE RECALIBRATED BAM FILES
    if (params.no_intervals) {
        cram_recalibrated = cram_applybqsr
    } else {
        cram_applybqsr.map{ meta, cram ->
            meta.id = meta.sample
            [meta, cram]
        }.groupTuple(size: num_intervals).set{cram_recalibrated_interval}

        SAMTOOLS_MERGE_CRAM(cram_recalibrated_interval, fasta)
        cram_recalibrated = SAMTOOLS_MERGE_CRAM.out.cram
        ch_versions = ch_versions.mix(SAMTOOLS_MERGE_CRAM.out.versions)
    }

    INDEX_RECALIBRATE(cram_recalibrated)
    cram_recalibrated_index = INDEX_RECALIBRATE.out.cram_crai
    ch_versions = ch_versions.mix(INDEX_RECALIBRATE.out.versions)

    qualimap_bamqc = Channel.empty()
    samtools_stats = Channel.empty()

    if (!skip_bamqc) {

        if(!params.wes || params.no_intervals) intervals_combined_bed_gz_tbi = [] //TODO: intervals also with WGS data? Probably need a parameter if WGS for deepvariant tool, that would allow to check here too

        QUALIMAP_BAMQC_CRAM(cram_recalibrated_index, intervals_combined_bed_gz_tbi, fasta, fasta_fai)
        qualimap_bamqc = QUALIMAP_BAMQC_CRAM.out.results
        ch_versions = ch_versions.mix(QUALIMAP_BAMQC_CRAM.out.versions)
    }

    if (!skip_samtools) {
        SAMTOOLS_STATS(cram_recalibrated_index, fasta)
        samtools_stats = SAMTOOLS_STATS.out.stats
        ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)
    }
    cram_reports = samtools_stats.mix(qualimap_bamqc)


    emit:
        cram = cram_recalibrated_index
        qc  = cram_reports
        versions = ch_versions
}
