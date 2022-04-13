//
// RECALIBRATE SPARK
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_APPLYBQSR_SPARK as APPLYBQSR_SPARK } from '../../../../modules/nf-core/modules/gatk4/applybqsrspark/main'
include { MERGE_INDEX_CRAM                         } from '../../merge_index_cram'

workflow RECALIBRATE_SPARK {
    take:
        cram          // channel: [mandatory] cram
        dict          // channel: [mandatory] dict
        fasta         // channel: [mandatory] fasta
        fasta_fai     // channel: [mandatory] fasta_fai
        intervals     // channel: [mandatory] intervals
        num_intervals //   value: [mandatory] number of intervals

    main:
    ch_versions = Channel.empty()

    cram_intervals = cram.combine(intervals)
        .map{ meta, cram, crai, recal, intervals ->
            new_meta = meta.clone()
            new_meta.id = num_intervals == 1 ? meta.sample : meta.sample + "_" + intervals.baseName
            [new_meta, cram, crai, recal, intervals]
        }

    // Run Applybqsr spark
    APPLYBQSR_SPARK(cram_intervals, fasta, fasta_fai, dict)

    // STEP 4.5: MERGING AND INDEXING THE RECALIBRATED BAM FILES
    MERGE_INDEX_CRAM(APPLYBQSR_SPARK.out.cram, fasta, num_intervals)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(APPLYBQSR_SPARK.out.versions)
    ch_versions = ch_versions.mix(MERGE_INDEX_CRAM.out.versions)

    emit:
        cram     = MERGE_INDEX_CRAM.out.cram_crai

        versions = ch_versions // channel: [ versions.yml ]
}
