//
// RECALIBRATE SPARK
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4SPARK_APPLYBQSR      } from '../../../modules/nf-core/gatk4spark/applybqsr'
include { CRAM_MERGE_INDEX_SAMTOOLS } from '../cram_merge_index_samtools'

workflow BAM_APPLYBQSR_SPARK {
    take:
    cram      // channel: [mandatory] [ meta, cram, crai, recal ]
    dict      // channel: [mandatory] [ dict ]
    fasta     // channel: [mandatory] [ fasta ]
    fasta_fai // channel: [mandatory] [ fasta_fai ]
    intervals // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals

    main:
    versions = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    // Move num_intervals to meta map
    cram_intervals = cram
        .combine(intervals)
        .map { meta, cram_, crai, recal, intervals_, num_intervals -> [meta + [num_intervals: num_intervals], cram_, crai, recal, intervals_] }

    // RUN APPLYBQSR SPARK
    GATK4SPARK_APPLYBQSR(
        cram_intervals,
        fasta.map { _meta, fasta_ -> [fasta_] },
        fasta_fai.map { _meta, fasta_fai_ -> [fasta_fai_] },
        dict.map { _meta, dict_ -> [dict_] },
    )

    // Gather the recalibrated cram files
    cram_to_merge = GATK4SPARK_APPLYBQSR.out.cram.map { meta, cram_ -> [groupKey(meta, meta.num_intervals), cram_] }.groupTuple()

    // Merge and index the recalibrated cram files
    CRAM_MERGE_INDEX_SAMTOOLS(
        cram_to_merge,
        fasta.map { _meta, fasta_ -> [fasta_] },
        fasta_fai.map { _meta, fasta_fai_ -> [fasta_fai_] },
    )

    // Remove no longer necessary field: num_intervals
    cram_recal = CRAM_MERGE_INDEX_SAMTOOLS.out.cram_crai.map { meta, cram_, crai -> [meta - meta.subMap('num_intervals'), cram_, crai] }

    // Gather versions of all tools used
    versions = versions.mix(GATK4SPARK_APPLYBQSR.out.versions)
    versions = versions.mix(CRAM_MERGE_INDEX_SAMTOOLS.out.versions)

    emit:
    cram     = cram_recal // channel: [ meta, cram, crai ]
    versions // channel: [ versions.yml ]
}
