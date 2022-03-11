//
// RECALIBRATE SPARK
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_APPLYBQSR_SPARK as APPLYBQSR_SPARK  } from '../../../../modules/local/gatk4/applybqsrspark/main'
include { SAMTOOLS_INDEX as INDEX_RECALIBRATE       } from '../../../../modules/local/samtools/index/main'
include { SAMTOOLS_MERGE_CRAM                       } from '../../../../modules/local/samtools/mergecram/main'

workflow RECALIBRATE_SPARK {
    take:
        cram           // channel: [mandatory] cram
        dict           // channel: [mandatory] dict
        fasta          // channel: [mandatory] fasta
        fasta_fai      // channel: [mandatory] fasta_fai
        intervals      // channel: [mandatory] intervals
        num_intervals
        no_intervals
        intervals_combined_bed_gz_tbi

    main:
    ch_versions = Channel.empty()
    qc_reports  = Channel.empty()

    cram_intervals = cram.combine(intervals)
        .map{ meta, cram, crai, recal, intervals ->
            new_meta = meta.clone()
            new_meta.id = intervals.baseName != "no_intervals" ? meta.sample + "_" + intervals.baseName : meta.sample
            [new_meta, cram, crai, recal, intervals]
        }

    // Run Applybqsr spark
    APPLYBQSR_SPARK(cram_intervals, fasta, fasta_fai, dict)

    cram_recalibrated_no_intervals = APPLYBQSR_SPARK.out.cram

    // Empty the no intervals cram channel if we have intervals
    if (!no_intervals) cram_recalibrated_no_intervals = Channel.empty()

    // STEP 4.5: MERGING AND INDEXING THE RECALIBRATED BAM FILES
    cram_recalibrated_interval = APPLYBQSR_SPARK.out.cram
        .map{ meta, cram ->
            meta.id = meta.sample
            [meta, cram]
        }.groupTuple(size: num_intervals)

    // Only when we have intervals
    SAMTOOLS_MERGE_CRAM(cram_recalibrated_interval, fasta)

    INDEX_RECALIBRATE(cram_recalibrated_no_intervals.mix(SAMTOOLS_MERGE_CRAM.out.cram))

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(APPLYBQSR_SPARK.out.versions)
    ch_versions = ch_versions.mix(INDEX_RECALIBRATE.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE_CRAM.out.versions)

    emit:
        cram     = INDEX_RECALIBRATE.out.cram_crai
        qc       = qc_reports

        versions = ch_versions // channel: [ versions.yml ]
}
