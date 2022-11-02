//
// RECALIBRATE
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_APPLYBQSR           } from '../../../modules/nf-core/gatk4/applybqsr/main'
include { CRAM_MERGE_INDEX_SAMTOOLS } from '../cram_merge_index_samtools/main'

workflow BAM_APPLYBQSR {
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
                data_type:      meta.data_type,
                id:             meta.sample,
                num_intervals:  num_intervals,
                patient:        meta.patient,
                sample:         meta.sample,
                sex:            meta.sex,
                status:         meta.status,
            ],
            cram, crai, recal, intervals_new]
        }

    // Run Applybqsr
    GATK4_APPLYBQSR(cram_intervals, fasta, fasta_fai, dict)

    // STEP 4.5: MERGING AND INDEXING THE RECALIBRATED CRAM FILES
    CRAM_MERGE_INDEX_SAMTOOLS(GATK4_APPLYBQSR.out.cram, fasta, fasta_fai)

    ch_cram_recal_out = CRAM_MERGE_INDEX_SAMTOOLS.out.cram_crai.map{ meta, cram, crai ->
                            // remove no longer necessary fields to make sure joining can be done correctly: num_intervals
                            [[
                                data_type:  meta.data_type,
                                id:         meta.id,
                                patient:    meta.patient,
                                sample:     meta.sample,
                                sex:        meta.sex,
                                status:     meta.status,
                            ],
                            cram, crai]
                        }

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(GATK4_APPLYBQSR.out.versions)
    ch_versions = ch_versions.mix(CRAM_MERGE_INDEX_SAMTOOLS.out.versions)

    emit:
        cram     = ch_cram_recal_out
        versions = ch_versions // channel: [ versions.yml ]
}
