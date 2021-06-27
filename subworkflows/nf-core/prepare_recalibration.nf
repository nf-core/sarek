/*
========================================================================================
    PREPARE RECALIBRATION
========================================================================================
*/

params.baserecalibrator_options  = [:]
params.gatherbqsrreports_options = [:]
params.baserecalibrator_spark_options  = [:]

include { GATK4_BASERECALIBRATOR  as BASERECALIBRATOR }              from '../../modules/nf-core/software/gatk4/baserecalibrator/main'  addParams(options: params.baserecalibrator_options)
include { GATK4_BASERECALIBRATOR_SPARK  as BASERECALIBRATOR_SPARK }  from '../../modules/nf-core/software/gatk4/baserecalibratorspark/main'  addParams(options: params.baserecalibrator_spark_options)
include { GATK4_GATHERBQSRREPORTS as GATHERBQSRREPORTS }             from '../../modules/nf-core/software/gatk4/gatherbqsrreports/main' addParams(options: params.gatherbqsrreports_options)

workflow PREPARE_RECALIBRATION {
    take:
        cram_markduplicates // channel: [mandatory] cram_markduplicates
        use_gatk_spark      //   value: [mandatory] use gatk spark
        dict                // channel: [mandatory] dict
        fai                 // channel: [mandatory] fai
        fasta               // channel: [mandatory] fasta
        intervals           // channel: [mandatory] intervals
        known_sites         // channel: [optional]  known_sites
        known_sites_tbi     // channel: [optional]  known_sites_tbi
        no_intervals        //   value: [mandatory] no_intervals
        known_indels
        dbsnp

    main:
    intervals.dump(tag:'intervals2')
    cram_markduplicates.combine(intervals).map{ meta, cram, crai, intervals ->
        new_meta = meta.clone()
        new_meta.id = meta.sample + "_" + intervals.baseName
        [new_meta, cram, crai, intervals]
    }.set{cram_markduplicates_intervals}

    cram_markduplicates_intervals.dump(tag:'bqsrinput')

    if(use_gatk_spark){
        BASERECALIBRATOR_SPARK(cram_markduplicates_intervals, fasta, fai, dict, known_sites, known_sites_tbi)
        table_baserecalibrator = BASERECALIBRATOR_SPARK.out.table
    }else{
        BASERECALIBRATOR(cram_markduplicates_intervals, fasta, fai, dict, known_sites_tbi, known_indels, dbsnp)
        table_baserecalibrator = BASERECALIBRATOR.out.table
    }

    // STEP 3.5: MERGING RECALIBRATION TABLES
    if (no_intervals) {
        table_baserecalibrator.map { meta, table ->
            meta.id = meta.sample
            [meta, table]
        }.set{table_bqsr}
    } else {
        table_baserecalibrator.map{ meta, table ->
            meta.id = meta.sample
            [meta, table]
        }.groupTuple().set{recaltable}


        GATHERBQSRREPORTS(recaltable)
        table_bqsr = GATHERBQSRREPORTS.out.table

    }

    emit:
        table_bqsr = table_bqsr
}
