//
// RECALIBRATE SPARK
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BAM_MERGE_INDEX_SAMTOOLS  } from '../bam_merge_index_samtools'
include { CRAM_MERGE_INDEX_SAMTOOLS } from '../cram_merge_index_samtools'
include { GATK4SPARK_APPLYBQSR      } from '../../../modules/nf-core/gatk4spark/applybqsr'

workflow BAM_APPLYBQSR_SPARK {
    take:
    cram      // channel: [mandatory] [ meta, cram, crai, recal ]
    dict      // channel: [mandatory] [ dict ]
    fasta     // channel: [mandatory] [ fasta ]
    fasta_fai // channel: [mandatory] [ fasta_fai ]
    intervals // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals

    main:
    versions = Channel.empty()
    bam_applybqsr_single = Channel.empty()
    bam_to_merge = Channel.empty()

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

    // FOR BAMs
    if (params.save_output_as_bam) {

        bam_applybqsr_out = GATK4SPARK_APPLYBQSR.out.bam
            .join(GATK4SPARK_APPLYBQSR.out.bai, failOnDuplicate: true, failOnMismatch: true)
            .branch {
                single: it[0].num_intervals == 1
                multiple: it[0].num_intervals > 1
            }

        bam_applybqsr_single = bam_applybqsr_out.single

        // For multiple intervals, gather and merge the recalibrated cram files
        bam_to_merge = bam_applybqsr_out.multiple
            .map { meta, bam_, _bai -> [groupKey(meta, meta.num_intervals), bam_] }
            .groupTuple()
    }

    // Merge and index the recalibrated cram files
    BAM_MERGE_INDEX_SAMTOOLS(bam_to_merge)

    // Combine single and merged multiple bam and index files, removing num_intervals field
    bam_recal = bam_applybqsr_single
        .mix(BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai)
        .map { meta, bam, bai -> [meta - meta.subMap('num_intervals'), bam, bai] }

    // FOR CRAMs

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
    versions = versions.mix(BAM_MERGE_INDEX_SAMTOOLS.out.versions)
    versions = versions.mix(CRAM_MERGE_INDEX_SAMTOOLS.out.versions)
    versions = versions.mix(GATK4SPARK_APPLYBQSR.out.versions)

    emit:
    bam      = bam_recal // channel: [ meta, bam, bai ]
    cram     = cram_recal // channel: [ meta, cram, crai ]
    versions // channel: [ versions.yml ]
}
