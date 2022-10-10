//
// RECALIBRATE SPARK
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_APPLYBQSR_SPARK as APPLYBQSR_SPARK } from '../../../../modules/nf-core/gatk4/applybqsrspark/main'
include { MERGE_INDEX_CRAM                         } from '../../merge_index_cram'

workflow RECALIBRATE_SPARK {
    take:
        cram          // channel: [mandatory] meta, cram, crai, recal
        dict          // channel: [mandatory] dict
        fasta         // channel: [mandatory] fasta
        fasta_fai     // channel: [mandatory] fasta_fai
        intervals     // channel: [mandatory] intervals, num_intervals

    main:
    ch_versions = Channel.empty()

    cram_intervals = cram.combine(intervals)
        .map{ meta, cram, crai, recal, intervals, num_intervals ->

            //If no interval file provided (0) then add empty list
            intervals_new = num_intervals == 0 ? [] : intervals

            [[
                id:             meta.sample,
                data_type:      meta.data_type,
                num_intervals:  num_intervals,
                patient:        meta.patient,
                sample:         meta.sample,
                sex:            meta.sex,
                status:         meta.status,
            ],
            cram, crai, recal, intervals_new]
        }

    // Run Applybqsr spark
    APPLYBQSR_SPARK(cram_intervals, fasta, fasta_fai, dict)

    // STEP 4.5: MERGING AND INDEXING THE RECALIBRATED BAM FILES
    MERGE_INDEX_CRAM(APPLYBQSR_SPARK.out.cram, fasta, fasta_fai)

    ch_cram_recal_out = MERGE_INDEX_CRAM.out.cram_crai.map{ meta, cram, crai ->
                             // remove no longer necessary fields to make sure joining can be done correctly: num_intervals
                            [[
                                id:         meta.id,
                                data_type:  meta.data_type,
                                patient:    meta.patient,
                                sample:     meta.sample,
                                sex:        meta.sex,
                                status:     meta.status,
                            ],
                            cram, crai]
                        }

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(APPLYBQSR_SPARK.out.versions)
    ch_versions = ch_versions.mix(MERGE_INDEX_CRAM.out.versions)

    emit:
        cram     = ch_cram_recal_out
        versions = ch_versions // channel: [ versions.yml ]
}
