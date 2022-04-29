//
// RECALIBRATE
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_APPLYBQSR as APPLYBQSR } from '../../../../modules/nf-core/modules/gatk4/applybqsr/main'
include { MERGE_INDEX_CRAM             } from '../../merge_index_cram'

workflow RECALIBRATE {
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
            new_meta = meta.clone()

            // If either no scatter/gather is done, i.e. no interval (0) or one interval (1), then don't rename samples
            new_meta.id = num_intervals <= 1 ? meta.sample : meta.sample + "_" + intervals.baseName
            new_meta.num_intervals = num_intervals

            //If no interval file provided (0) then add empty list
            intervals_new = num_intervals == 0 ? [] : intervals

            [new_meta, cram, crai, recal, intervals_new]
        }

    // Run Applybqsr
    APPLYBQSR(cram_intervals, fasta, fasta_fai, dict)

    // STEP 4.5: MERGING AND INDEXING THE RECALIBRATED CRAM FILES
    MERGE_INDEX_CRAM(APPLYBQSR.out.cram, fasta)

    ch_cram_recal_out = MERGE_INDEX_CRAM.out.cram_crai.map{ meta, cram, crai ->
                            new_meta = meta.clone()

                            // remove no longer necessary fields to make sure joining can be done correctly
                            new_meta.remove('num_intervals')

                            [new_meta, cram, crai]
                        }

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(APPLYBQSR.out.versions)
    ch_versions = ch_versions.mix(MERGE_INDEX_CRAM.out.versions)

    emit:
        cram     = ch_cram_recal_out
        versions = ch_versions // channel: [ versions.yml ]
}
