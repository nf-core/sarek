include { GATK4_APPLYBQSR } from '../../../modules/nf-core/gatk4/applybqsr'
include { SAMTOOLS_MERGE  } from '../../../modules/nf-core/samtools/merge'

workflow BAM_APPLYBQSR {
    take:
    input // channel: [mandatory] [ meta, reads, index, recal ]
    references // channel: [mandatory] [ meta, fasta, fai, dict ]
    intervals // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals
    output_format // tuple: [mandatory] [ file_type, index_format ] — ['bam','bai'], ['bam','csi'], or ['cram','crai']

    main:
    reads_interval = input
        .combine(intervals)
        .map { meta, reads_, index, recal, intervals_, num_intervals ->
            [meta + [num_intervals: num_intervals, interval: intervals_.name], reads_, index, recal, intervals_]
        }

    // RUN APPLYBQSR
    GATK4_APPLYBQSR(reads_interval, references, output_format)

    reads_applybqsr = GATK4_APPLYBQSR.out.reads
        .join(GATK4_APPLYBQSR.out.index, failOnMismatch: true)
        .branch { meta, _reads, _index ->
            single: meta.num_intervals <= 1
            multiple: meta.num_intervals > 1
        }

    // For multiple intervals, gather and merge the recalibrated files
    reads_to_merge = reads_applybqsr.multiple
        .map { meta, reads, index ->
            def merge_meta = meta - meta.subMap('interval')
            [groupKey(merge_meta, merge_meta.num_intervals), reads, index]
        }
        .groupTuple()

    SAMTOOLS_MERGE(reads_to_merge, references.map { meta, fasta, fai, _dict -> [meta, fasta, fai, []] }, output_format[1])

    // Unified output — bam or cram depending on what was produced
    recal_out = SAMTOOLS_MERGE.out.bam
        .mix(SAMTOOLS_MERGE.out.cram)
        .join(SAMTOOLS_MERGE.out.index, failOnMismatch: true)
        .mix(reads_applybqsr.single)
        .map { meta, reads, index -> [meta - meta.subMap('num_intervals', 'interval'), reads, index] }

    emit:
    recal_out = recal_out // channel: [ meta, reads, index ]
}
