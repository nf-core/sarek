//
// PREPARE RECALIBRATION with SPARK
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_BASERECALIBRATOR_SPARK as BASERECALIBRATOR_SPARK } from '../../../../modules/nf-core/modules/gatk4/baserecalibratorspark/main'
include { GATK4_GATHERBQSRREPORTS      as GATHERBQSRREPORTS      } from '../../../../modules/nf-core/modules/gatk4/gatherbqsrreports/main'

workflow PREPARE_RECALIBRATION_SPARK {
    take:
        cram            // channel: [mandatory] cram_markduplicates
        dict            // channel: [mandatory] dict
        fasta           // channel: [mandatory] fasta
        fasta_fai       // channel: [mandatory] fasta_fai
        intervals       // channel: [mandatory] intervals
        known_sites     // channel: [optional]  known_sites
        known_sites_tbi // channel: [optional]  known_sites_tbi

    main:
    ch_versions = Channel.empty()

    cram_intervals = cram.combine(intervals)
        .map{ meta, cram, crai, intervals, num_intervals ->
            new_meta = meta.clone()

            // If either no scatter/gather is done, i.e. no interval (0) or one interval (1), then don't rename samples
            new_meta.id = num_intervals <= 1 ? meta.sample : meta.sample + "_" + intervals.baseName
            new_meta.num_intervals = num_intervals

            //If no interval file provided (0) then add empty list
            intervals_new = num_intervals == 0 ? [] : intervals

            [new_meta, cram, crai, intervals_new]
        }

    // Run Baserecalibrator spark
    BASERECALIBRATOR_SPARK(cram_intervals, fasta, fasta_fai, dict, known_sites, known_sites_tbi)

    // Figuring out if there is one or more table(s) from the same sample
    table_to_merge = BASERECALIBRATOR_SPARK.out.table
        .map{ meta, table ->
                meta.id = meta.sample

                def groupKey = groupKey(meta, meta.num_intervals)
                [meta, table]
        }.groupTuple()
    .branch{
        single:   it[1].size() == 1
        multiple: it[1].size() > 1
    }

    // STEP 3.5: MERGING RECALIBRATION TABLES

    // Merge the tables only when we have intervals
    GATHERBQSRREPORTS(table_to_merge.multiple)
    table_bqsr = table_to_merge.single.mix(GATHERBQSRREPORTS.out.table)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(BASERECALIBRATOR_SPARK.out.versions)
    ch_versions = ch_versions.mix(GATHERBQSRREPORTS.out.versions)

    emit:
        table_bqsr = table_bqsr
        versions   = ch_versions // channel: [versions.yml]
}
