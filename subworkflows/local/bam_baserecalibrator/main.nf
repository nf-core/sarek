//
// PREPARE RECALIBRATION
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_BASERECALIBRATOR  } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_GATHERBQSRREPORTS } from '../../../modules/nf-core/gatk4/gatherbqsrreports/main'

workflow BAM_BASERECALIBRATOR {
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

            [[
                data_type:      meta.data_type,
                id:             meta.sample,
                num_intervals:  num_intervals,
                patient:        meta.patient,
                sample:         meta.sample,
                sex:            meta.sex,
                status:         meta.status,
            ],
            cram, crai, intervals_new]
        }

    // Run Baserecalibrator
    GATK4_BASERECALIBRATOR(cram_intervals, fasta, fasta_fai, dict, known_sites, known_sites_tbi)

    // Figuring out if there is one or more table(s) from the same sample
    table_to_merge = GATK4_BASERECALIBRATOR.out.table
        .map{ meta, table ->

                new_meta = [
                                data_type:      meta.data_type,
                                id:             meta.sample,
                                num_intervals:  meta.num_intervals,
                                patient:        meta.patient,
                                sample:         meta.sample,
                                sex:            meta.sex,
                                status:         meta.status,
                            ]

                [groupKey(new_meta, meta.num_intervals), table]
        }.groupTuple()
    .branch{
        //Warning: size() calculates file size not list length here, so use num_intervals instead
        single:   it[0].num_intervals <= 1
        multiple: it[0].num_intervals > 1
    }

    // STEP 3.5: MERGING RECALIBRATION TABLES

    // Merge the tables only when we have intervals
    GATK4_GATHERBQSRREPORTS(table_to_merge.multiple)
    table_bqsr = table_to_merge.single.map{meta, table -> [meta, table[0]]}.mix(GATK4_GATHERBQSRREPORTS.out.table)
                                        .map{ meta, table ->
                                            // remove no longer necessary fields to make sure joining can be done correctly: num_intervals
                                            [[
                                                data_type:  meta.data_type,
                                                id:         meta.sample,
                                                patient:    meta.patient,
                                                sample:     meta.sample,
                                                sex:        meta.sex,
                                                status:     meta.status,
                                            ],
                                            table]
                                        }

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions)
    ch_versions = ch_versions.mix(GATK4_GATHERBQSRREPORTS.out.versions)

    emit:
        table_bqsr = table_bqsr
        versions   = ch_versions // channel: [versions.yml]
}
