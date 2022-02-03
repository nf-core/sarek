//
// PREPARE RECALIBRATION
//

include { GATK4_BASERECALIBRATOR as BASERECALIBRATOR             } from '../../modules/nf-core/modules/gatk4/baserecalibrator/main'
include { GATK4_BASERECALIBRATOR_SPARK as BASERECALIBRATOR_SPARK } from '../../modules/local/gatk4/baserecalibratorspark/main'
include { GATK4_GATHERBQSRREPORTS as GATHERBQSRREPORTS           } from '../../modules/nf-core/modules/gatk4/gatherbqsrreports/main'

workflow PREPARE_RECALIBRATION {
    take:
        cram_markduplicates // channel: [mandatory] cram_markduplicates
        use_gatk_spark      //   value: [mandatory] use gatk spark
        dict                // channel: [mandatory] dict
        fasta               // channel: [mandatory] fasta
        fasta_fai           // channel: [mandatory] fasta_fai
        intervals           // channel: [mandatory] intervals
        num_intervals
        known_sites         // channel: [optional]  known_sites
        known_sites_tbi     // channel: [optional]  known_sites_tbi
        no_intervals        //   value: [mandatory] no_intervals

    main:

    ch_versions = Channel.empty()

    cram_markduplicates.combine(intervals)
        .map{ meta, cram, crai, intervals ->
            new_meta = meta.clone()
            new_meta.id = intervals.baseName != "no_intervals" ? meta.sample + "_" + intervals.baseName : meta.sample
            [new_meta, cram, crai, intervals]
        }.set{cram_markduplicates_intervals}

    if (use_gatk_spark) {
        BASERECALIBRATOR_SPARK(cram_markduplicates_intervals, fasta, fasta_fai, dict, known_sites, known_sites_tbi)
        table_baserecalibrator = BASERECALIBRATOR_SPARK.out.table
        ch_versions = ch_versions.mix(BASERECALIBRATOR_SPARK.out.versions)

    } else {
        BASERECALIBRATOR(cram_markduplicates_intervals, fasta, fasta_fai, dict, known_sites, known_sites_tbi)
        table_baserecalibrator = BASERECALIBRATOR.out.table
        ch_versions = ch_versions.mix(BASERECALIBRATOR.out.versions)
    }

    //STEP 3.5: MERGING RECALIBRATION TABLES
    if (no_intervals) {
        table_baserecalibrator.map { meta, table ->
            meta.id = meta.sample
            [meta, table]
        }.set{table_bqsr}
    } else {
        table_baserecalibrator.map{ meta, table ->
            meta.id = meta.sample
            [meta, table]
        }.groupTuple(size: num_intervals).set{recaltable}

        GATHERBQSRREPORTS(recaltable)
        table_bqsr = GATHERBQSRREPORTS.out.table
        ch_versions = ch_versions.mix(GATHERBQSRREPORTS.out.versions)
    }

    emit:
        table_bqsr = table_bqsr
        versions   = ch_versions // channel: [versions.yml]
}
