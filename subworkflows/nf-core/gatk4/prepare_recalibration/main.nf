//
// PREPARE RECALIBRATION
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_BASERECALIBRATOR       as BASERECALIBRATOR       } from '../../../../modules/nf-core/modules/gatk4/baserecalibrator/main'
include { GATK4_BASERECALIBRATOR_SPARK as BASERECALIBRATOR_SPARK } from '../../../../modules/local/gatk4/baserecalibratorspark/main'
include { GATK4_GATHERBQSRREPORTS      as GATHERBQSRREPORTS      } from '../../../../modules/nf-core/modules/gatk4/gatherbqsrreports/main'

workflow PREPARE_RECALIBRATION {
    take:
        cram_markduplicates // channel: [mandatory] cram_markduplicates
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

    cram_markduplicates_intervals = cram_markduplicates.combine(intervals)
        .map{ meta, cram, crai, intervals ->
            new_meta = meta.clone()
            new_meta.id = intervals.baseName != "no_intervals" ? meta.sample + "_" + intervals.baseName : meta.sample
            [new_meta, cram, crai, intervals]
        }

    // Run Baserecalibrator or Baserecalibrator spark
    BASERECALIBRATOR(cram_markduplicates_intervals, fasta, fasta_fai, dict, known_sites, known_sites_tbi)
    BASERECALIBRATOR_SPARK(cram_markduplicates_intervals, fasta, fasta_fai, dict, known_sites, known_sites_tbi)

    table_baserecalibrator = BASERECALIBRATOR.out.table.mix(BASERECALIBRATOR_SPARK.out.table)
        .map{ meta, table ->
                meta.id = meta.sample
                [meta, table]
            }

    // Only one of the two channels will be used
    table_no_intervals = table_baserecalibrator
    table_intervals    = table_baserecalibrator.groupTuple(size: num_intervals)

    // STEP 3.5: MERGING RECALIBRATION TABLES
    // Empty the no intervals table channel if we have intervals
    if (!no_intervals) table_no_intervals = Channel.empty()

    // Merge the tables only when we have intervals
    GATHERBQSRREPORTS(table_intervals)
    table_bqsr = table_no_intervals.mix(GATHERBQSRREPORTS.out.table)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(BASERECALIBRATOR.out.versions)
    ch_versions = ch_versions.mix(BASERECALIBRATOR_SPARK.out.versions)
    ch_versions = ch_versions.mix(GATHERBQSRREPORTS.out.versions)

    emit:
        table_bqsr = table_bqsr
        versions   = ch_versions // channel: [versions.yml]
}
