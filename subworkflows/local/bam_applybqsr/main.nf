//
// RECALIBRATE
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BAM_MERGE_INDEX_SAMTOOLS  } from '../bam_merge_index_samtools'
include { CRAM_MERGE_INDEX_SAMTOOLS } from '../cram_merge_index_samtools'
include { GATK4_APPLYBQSR           } from '../../../modules/nf-core/gatk4/applybqsr'

workflow BAM_APPLYBQSR {
    take:
    cram      // channel: [mandatory] [ meta, cram, crai, recal ]
    dict      // channel: [mandatory] [ meta, dict ]
    fasta     // channel: [mandatory] [ meta, fasta ]
    fasta_fai // channel: [mandatory] [ meta, fasta_fai ]
    intervals // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals

    main:
    versions = channel.empty()

    // Combine cram and intervals for spread and gather strategy
    // Move num_intervals to meta map
    cram_intervals = cram
        .combine(intervals)
        .map { meta, cram_, crai, recal, intervals_, num_intervals -> [meta + [num_intervals: num_intervals], cram_, crai, recal, intervals_] }

    // RUN APPLYBQSR
    GATK4_APPLYBQSR(
        cram_intervals,
        fasta.map { _meta, fasta_ -> [fasta_] },
        fasta_fai.map { _meta, fasta_fai_ -> [fasta_fai_] },
        dict.map { _meta, dict_ -> [dict_] },
    )

    // BAM path — populated when ext.suffix='bam', empty otherwise
    bam_to_merge = GATK4_APPLYBQSR.out.bam
        .map { meta, bam_ -> [groupKey(meta, meta.num_intervals), bam_] }
        .groupTuple()

    BAM_MERGE_INDEX_SAMTOOLS(bam_to_merge)

    // CRAM path — populated when ext.suffix='cram', empty otherwise
    cram_to_merge = GATK4_APPLYBQSR.out.cram.map { meta, cram_ -> [groupKey(meta, meta.num_intervals), cram_] }.groupTuple()

    CRAM_MERGE_INDEX_SAMTOOLS(
        cram_to_merge,
        fasta,
        fasta_fai,
    )

    // Mix — one is always empty
    recal_out = BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai
        .mix(CRAM_MERGE_INDEX_SAMTOOLS.out.cram_crai)
        .map { meta, file_, index -> [meta - meta.subMap('num_intervals'), file_, index] }

    // Gather versions of all tools used
    versions = versions.mix(BAM_MERGE_INDEX_SAMTOOLS.out.versions)
    versions = versions.mix(CRAM_MERGE_INDEX_SAMTOOLS.out.versions)
    versions = versions.mix(GATK4_APPLYBQSR.out.versions)

    emit:
    alignment = recal_out // channel: [ meta, file, index ] — BAM or CRAM
    versions              // channel: [ versions.yml ]
}
