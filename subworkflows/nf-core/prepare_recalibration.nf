/*
========================================================================================
    PREPARE RECALIBRATION
========================================================================================
*/

params.baserecalibrator_options  = [:]
params.gatherbqsrreports_options = [:]

include { GATK4_BASERECALIBRATOR  as BASERECALIBRATOR }  from '../../modules/nf-core/software/gatk4/baserecalibrator/main'  addParams(options: params.baserecalibrator_options)
include { GATK4_GATHERBQSRREPORTS as GATHERBQSRREPORTS } from '../../modules/nf-core/software/gatk4/gatherbqsrreports/main' addParams(options: params.gatherbqsrreports_options)

workflow PREPARE_RECALIBRATION {
    take:
        bam_markduplicates // channel: [mandatory] bam_markduplicates
        dict               // channel: [mandatory] dict
        fai                // channel: [mandatory] fai
        fasta              // channel: [mandatory] fasta
        intervals          // channel: [mandatory] intervals
        known_sites        // channel: [optional]  known_sites
        known_sites_tbi    // channel: [optional]  known_sites_tbi
        no_intervals       //   value: [mandatory] no_intervals

    main:

    bam_markduplicates.combine(intervals).map{ meta, bam, bai, intervals ->
        new_meta = meta.clone()
        new_meta.id = meta.sample + "_" + intervals.baseName
        [new_meta, bam, bai, intervals]
    }.set{bam_markduplicates_intervals}

    BASERECALIBRATOR(bam_markduplicates_intervals, fasta, fai, dict, known_sites, known_sites_tbi)

    // STEP 3.5: MERGING RECALIBRATION TABLES
    if (no_intervals) {
        BASERECALIBRATOR.out.table.map { meta, table ->
            meta.id = meta.sample
            [meta, table]
        }.set{table_bqsr}
    } else {
        BASERECALIBRATOR.out.table
            .map{ meta, table ->
            meta.id = meta.sample
            [meta, table]
        }.groupTuple().set{recaltable}

        GATHERBQSRREPORTS(recaltable)
        table_bqsr = GATHERBQSRREPORTS.out.table
    }

    emit:
        table_bqsr = table_bqsr
}
