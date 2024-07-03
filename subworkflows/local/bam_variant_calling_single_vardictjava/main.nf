//
// Vardictjava germline calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { VARDICTJAVA                                  } from '../../../modules/nf-core/vardictjava/main'
include { GATK4_MERGEVCFS  as MERGE_VARDICTJAVA        } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM              } from '../../../modules/nf-core/samtools/convert/main'

workflow BAM_VARIANT_CALLING_SINGLE_VARDICTJAVA {
    take:
    cram          // channel: [mandatory] [ meta, cram, crai ]
    dict          // channel: [optional]  [ meta, dict ]
    fasta         // channel: [mandatory] [ fasta ]
    fasta_fai     // channel: [mandatory] [ fasta_fai ]
    intervals     // channel: [mandatory] [ interval.bed.gz, interval.bed.gz.tbi, num_intervals ] or [ [], [], 0 ] if no intervals

    main:
    versions = Channel.empty()

    // Convert cram to bam
    cram
    .branch {meta, cram, crai ->
        bam:    cram.extension == "bam"
        cram:   cram.extension == "cram"}
        .set{ch_bam_from_cram}

    CRAM_TO_BAM(
        ch_bam_from_cram.cram,
        fasta,
        fasta_fai
    )

    // Combine converted bam, bai and intervals
    ch_bam_from_cram.bam
        .mix(CRAM_TO_BAM.out.alignment_index)
        .combine(intervals)
        .map{meta, bam, bai, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], bam, bai, intervals ]}
        .set{ ch_vardict_input}

    VARDICTJAVA(
        ch_vardict_input,
        fasta.map{fasta -> [[id:fasta.baseName], fasta]},
        fasta_fai.map{fasta_fai -> [[id:fasta_fai.baseName], fasta_fai]}
    )

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf = VARDICTJAVA.out.vcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    vcf_to_merge        = vcf.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()

    MERGE_VARDICTJAVA(vcf_to_merge, dict)

    // Mix intervals and no_intervals channels together
    vcf = Channel.empty().mix(MERGE_VARDICTJAVA.out.vcf, vcf.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'vardictjava' ], vcf ] }

    versions = versions.mix(VARDICTJAVA.out.versions)
    versions = versions.mix(MERGE_VARDICTJAVA.out.versions)

    emit:
    vcf

    versions
}
