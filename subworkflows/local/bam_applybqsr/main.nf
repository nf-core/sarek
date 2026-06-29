//
// RECALIBRATE
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAMTOOLS_MERGE as MERGE_BAM  } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_MERGE as MERGE_CRAM } from '../../../modules/nf-core/samtools/merge/main'
include { GATK4_APPLYBQSR as GATK4_APPLYBQSR_CRAM } from '../../../modules/nf-core/gatk4/applybqsr'
include { GATK4_APPLYBQSR as GATK4_APPLYBQSR_BAM } from '../../../modules/nf-core/gatk4/applybqsr'

workflow BAM_APPLYBQSR {
    take:
    bam       // channel: [mandatory] [ meta, bam, bai, recal ]
    cram      // channel: [mandatory] [ meta, cram, crai, recal ]
    fasta     // channel: [mandatory] [ meta, fasta ]
    fasta_fai // channel: [mandatory] [ meta, fasta_fai ]
    dict      // channel: [mandatory] [ dict ]
    intervals // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals

    main:
    bam_recal = Channel.empty()
    cram_recal = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    // Move num_intervals to meta map
    cram_intervals = cram
        .combine(intervals)
        .map { meta, cram_, crai, recal, intervals_, num_intervals -> [meta + [num_intervals: num_intervals], cram_, crai, recal, intervals_] }

    bam_intervals = bam
        .combine(intervals)
        .map { meta, bam_, bai, recal, intervals_, num_intervals -> [meta + [num_intervals: num_intervals], bam_, bai, recal, intervals_] }

    // RUN APPLYBQSR 
    GATK4_APPLYBQSR_CRAM(
        cram_intervals,
        fasta.map { _meta, fasta_file -> [fasta_file] },
        fasta_fai.map { _meta, fasta_file, fai_file -> [fasta_file, fai_file] },
        dict.map { _meta, dict_file -> [dict_file] },
    )

    // RUN APPLYBQSR 
    GATK4_APPLYBQSR_BAM(
        bam_intervals,
        fasta.map { _meta, fasta_file -> [fasta_file] },
        fasta_fai.map { _meta, fasta_file, fai_file -> [fasta_file, fai_file] },
        dict.map { _meta, dict_file -> [dict_file] },
    )

    // FOR BAMs
    if (params.save_output_as_bam) {

        bam_applybqsr_out = GATK4_APPLYBQSR_BAM.out.bam
            .join(GATK4_APPLYBQSR_BAM.out.bai, failOnDuplicate: true, failOnMismatch: true)
            .branch {
                single: it[0].num_intervals == 1
                multiple: it[0].num_intervals > 1
            }

        bam_applybqsr_single = bam_applybqsr_out.single

        // For multiple intervals, gather and merge the recalibrated cram files
        bam_to_merge = bam_applybqsr_out.multiple
            .map { meta, bam_, _bai -> [groupKey(meta, meta.num_intervals), bam_] }
            .groupTuple()

        MERGE_BAM(
                bam_to_merge,
                [ [ id:'null' ], [] ],
                [ [ id:'null' ], [] ]
            )
        
        bam_all = MERGE_BAM.out.bam.mix(bam_applybqsr_single)

    // Combine single and merged multiple bam and index files, removing num_intervals field
        bam_recal = bam_all
            .join(MERGE_BAM.out.bam_bai, failOnDuplicate: true, failOnMismatch: true)
            .map { meta, bam, bai -> [meta - meta.subMap('num_intervals'), bam, bai] }

    // FOR CRAMs

        // Branch CRAM output by number of intervals
    cram_applybqsr_out = GATK4_APPLYBQSR_CRAM.out.cram
        .join(GATK4_APPLYBQSR_CRAM.out.crai, failOnDuplicate: true, failOnMismatch: true)
        .branch {
            single: it[0].num_intervals == 1
            multiple: it[0].num_intervals > 1
        }
    
    cram_applybqsr_single = cram_applybqsr_out.single

    // Gather the recalibrated cram files
    cram_to_merge = cram_applybqsr_out.multiple
            .map { meta, cram_, _crai -> [groupKey(meta, meta.num_intervals), cram_] }
            .groupTuple()

    // Merge and index the recalibrated cram files
    MERGE_CRAM(
        cram_to_merge,
        [ [ id:'null' ], [] ],
        [ [ id:'null' ], [] ]
    )
    cram_all = MERGE_CRAM.out.bam

    // Remove no longer necessary field: num_intervals
    cram_recal = cram_all
        .join(MERGE_CRAM.out.bai, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, cram, crai -> [meta - meta.subMap('num_intervals'), cram, crai] }

    emit:
    bam      = bam_recal // channel: [ meta, bam, bai ]
    cram     = cram_recal // channel: [ meta, cram, crai ]
    }
}
