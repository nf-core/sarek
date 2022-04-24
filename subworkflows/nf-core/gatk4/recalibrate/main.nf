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
            new_meta.id = num_intervals == 1 ? meta.sample : meta.sample + "_" + intervals.baseName
            new_meta.num_intervals = num_intervals
            intervals_new = params.no_intervals ? [] : intervals
            [new_meta, cram, crai, recal, intervals_new]
        }

    // Run Applybqsr
    APPLYBQSR(cram_intervals, fasta, fasta_fai, dict)

    // STEP 4.5: MERGING AND INDEXING THE RECALIBRATED BAM FILES
    MERGE_INDEX_CRAM(APPLYBQSR.out.cram, fasta)

    ch_cram_recal_out = MERGE_INDEX_CRAM.out.cram_crai.map{ meta, cram, crai ->
                            new_meta = meta.clone()
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
