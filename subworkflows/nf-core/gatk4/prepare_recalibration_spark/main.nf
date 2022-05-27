//
// PREPARE RECALIBRATION with SPARK
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_BASERECALIBRATOR_SPARK as BASERECALIBRATOR_SPARK } from '../../../../modules/nf-core/modules/gatk4/baserecalibratorspark/main'
include { GATK4_GATHERBQSRREPORTS      as GATHERBQSRREPORTS      } from '../../../../modules/nf-core/modules/gatk4/gatherbqsrreports/main'

workflow PREPARE_RECALIBRATION_SPARK {
    take:
        cram            // channel: [mandatory] meta, cram_markduplicates, crai
        dict            // channel: [mandatory] dict
        fasta           // channel: [mandatory] fasta
        fasta_fai       // channel: [mandatory] fasta_fai
        intervals       // channel: [mandatory] intervals, num_intervals
        known_sites     // channel: [optional]  known_sites
        known_sites_tbi // channel: [optional]  known_sites_tbi

    main:
    ch_versions = Channel.empty()

    cram_intervals = cram.combine(intervals)
        .map{ meta, cram, crai, intervals, num_intervals ->

            //If no interval file provided (0) then add empty list
            intervals_new = num_intervals == 0 ? [] : intervals

            // If either no scatter/gather is done, i.e. no interval (0) or one interval (1), then don't rename samples
            new_id = num_intervals <= 1 ? meta.sample : meta.sample + "_" + intervals_new.baseName

            [[patient:meta.patient, sample:meta.sample, gender:meta.gender, status:meta.status, id:new_id, data_type:meta.data_type, num_intervals:num_intervals],
            cram, crai, intervals_new]
        }

    // Run Baserecalibrator spark
    BASERECALIBRATOR_SPARK(cram_intervals, fasta, fasta_fai, dict, known_sites, known_sites_tbi)

    // Figuring out if there is one or more table(s) from the same sample
    table_to_merge = BASERECALIBRATOR_SPARK.out.table
        .map{ meta, table ->

                new_meta = [patient:meta.patient, sample:meta.sample, gender:meta.gender, status:meta.status, id:meta.sample, data_type:meta.data_type, num_intervals:meta.num_intervals]

                [groupKey(new_meta, meta.num_intervals), table]
        }.groupTuple()
    .branch{
        //Warning: size() calculates file size not list length here, so use num_intervals instead
        single:   it[0].num_intervals <= 1
        multiple: it[0].num_intervals > 1
    }

    // STEP 3.5: MERGING RECALIBRATION TABLES

    // Merge the tables only when we have intervals
    GATHERBQSRREPORTS(table_to_merge.multiple)
    table_bqsr = table_to_merge.single.mix(GATHERBQSRREPORTS.out.table)
                                        .map{ meta, table ->
                                            // remove no longer necessary fields to make sure joining can be done correctly: num_intervals
                                            [[patient:meta.patient, sample:meta.sample, gender:meta.gender, status:meta.status, id:meta.sample, data_type:meta.data_type],
                                            table]
                                        }

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(BASERECALIBRATOR_SPARK.out.versions)
    ch_versions = ch_versions.mix(GATHERBQSRREPORTS.out.versions)

    emit:
        table_bqsr = table_bqsr
        versions   = ch_versions // channel: [versions.yml]
}
