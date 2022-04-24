//
// PREPARE RECALIBRATION
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_BASERECALIBRATOR  as BASERECALIBRATOR  } from '../../../../modules/nf-core/modules/gatk4/baserecalibrator/main'
include { GATK4_GATHERBQSRREPORTS as GATHERBQSRREPORTS } from '../../../../modules/nf-core/modules/gatk4/gatherbqsrreports/main'

workflow PREPARE_RECALIBRATION {
    take:
        cram            // channel: [mandatory] meta, cram_markduplicates, crai
        dict            // channel: [mandatory] dict
        fasta           // channel: [mandatory] fasta
        fasta_fai       // channel: [mandatory] fasta_fai
        intervals       // channel: [mandatory] intervals, num_intervals
        known_sites     // channel: [optional]  known_sites
        known_sites_tbi // channel: [optional]  known_sites_tbi
        //num_intervals   //   value: [mandatory] number of intervals

    main:
    ch_versions = Channel.empty()

    cram_intervals = cram.combine(intervals)
        .map{ meta, cram, crai, intervals, num_intervals ->
            new_meta = meta.clone()
            new_meta.id = num_intervals == 1 ? meta.sample : meta.sample + "_" + intervals.baseName
            new_meta.num_intervals = num_intervals
            intervals_new = params.no_intervals ? [] : intervals
            [new_meta, cram, crai, intervals_new]
        }

    // Run Baserecalibrator
    BASERECALIBRATOR(cram_intervals, fasta, fasta_fai, dict, known_sites, known_sites_tbi)

    // Figuring out if there is one or more table(s) from the same sample
    table_to_merge = BASERECALIBRATOR.out.table
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
                                        .map{ meta, table ->
                                            new_meta = meta.clone()
                                                // remove no longer necessary fields to make sure joining can be done correctly
                                            new_meta.remove('num_intervals')
                                            [new_meta, table]
                                        }

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(BASERECALIBRATOR.out.versions)
    ch_versions = ch_versions.mix(GATHERBQSRREPORTS.out.versions)

    emit:
        table_bqsr = table_bqsr
        versions   = ch_versions // channel: [versions.yml]
}
