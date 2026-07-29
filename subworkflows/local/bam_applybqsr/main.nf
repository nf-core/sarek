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
    cram                // channel: [mandatory] [ meta, cram, crai, recal ]
    dict                // channel: [mandatory] [ meta, dict ]
    fasta               // channel: [mandatory] [ meta, fasta ]
    fasta_fai           // channel: [mandatory] [ meta, fasta_fai ]
    intervals           // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals
    save_output_as_bam  // boolean: [mandatory] params.save_output_as_bam

    main:

    output_suffix = save_output_as_bam ? 'bam' : 'cram'

    // Combine cram and intervals for spread and gather strategy
    // Move num_intervals to meta map
    cram_intervals = cram
        .combine(intervals)
        .map { meta, cram_, crai, recal, intervals_, num_intervals -> [meta + [num_intervals: num_intervals], cram_, crai, recal, intervals_] }

    // fasta/fasta_fai/dict carry independent meta maps (computed vs user-supplied) so they
    // can't be joined by key — reshape into one tuple the way CRAM_SAMPLEQC does in sarek.nf
    fasta_fai_dict = fasta
        .combine(fasta_fai)
        .combine(dict)
        .map { meta_fasta, fasta_, _meta_fai, fai, _meta_dict, dict_ -> [meta_fasta, fasta_, fai, dict_] }
        .collect()

    // RUN APPLYBQSR
    GATK4_APPLYBQSR(
        cram_intervals,
        fasta_fai_dict,
        output_suffix,
    )

    // BAM path — populated when output_suffix == 'bam', empty otherwise
    bam_to_merge = GATK4_APPLYBQSR.out.bam
        .map { meta, bam_ -> [groupKey(meta, meta.num_intervals), bam_] }
        .groupTuple()

    BAM_MERGE_INDEX_SAMTOOLS(bam_to_merge)

    // CRAM path — populated when output_suffix == 'cram', empty otherwise
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

    emit:
    alignment = recal_out // channel: [ meta, file, index ] — BAM or CRAM
}
