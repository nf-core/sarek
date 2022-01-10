//
// RECALIBRATE
//

include { GATK4_APPLYBQSR as APPLYBQSR              } from '../../modules/nf-core/modules/gatk4/applybqsr/main'
include { GATK4_APPLYBQSR_SPARK as APPLYBQSR_SPARK  } from '../../modules/local/gatk4/applybqsrspark/main'
include { QUALIMAP_BAMQC_CRAM                       } from '../../modules/local/qualimap/bamqccram/main'
include { SAMTOOLS_INDEX as INDEX_RECALIBRATE       } from '../../modules/nf-core/modules/samtools/index/main'
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
        APPLYBQSR_SPARK(cram_intervals, fasta, fasta_fai, dict)
        cram_applybqsr = APPLYBQSR_SPARK.out.cram
    }else{
        APPLYBQSR(cram_intervals, fasta, fasta_fai, dict)
        cram_applybqsr = APPLYBQSR.out.cram
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

        INDEX_RECALIBRATE(cram_recalibrated)
        cram_recalibrated_index = INDEX_RECALIBRATE.out.cram_crai

        qualimap_bamqc = Channel.empty()
        samtools_stats = Channel.empty()

        if (!skip_bamqc) {
            QUALIMAP_BAMQC_CRAM(cram_recalibrated_index, target_bed, fasta, fasta_fai)
            qualimap_bamqc = QUALIMAP_BAMQC_CRAM.out.results
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
