include { GATK4_MERGEVCFS as MERGE_MANTA_DIPLOID      } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_MANTA_SMALL_INDELS } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_MANTA_SV           } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { MANTA_GERMLINE                              } from '../../../modules/nf-core/manta/germline/main'

// Seems to be the consensus on upstream modules implementation too
workflow BAM_VARIANT_CALLING_GERMLINE_MANTA {
    take:
    cram          // channel: [mandatory] [ meta, cram, crai ]
    dict          // channel: [optional]  [ meta, dict ]
    fasta         // channel: [mandatory] [ fasta ]
    fasta_fai     // channel: [mandatory] [ fasta_fai ]
    intervals     // channel: [mandatory] [ interval.bed.gz, interval.bed.gz.tbi, num_intervals ] or [ [], [], 0 ] if no intervals

    main:
    versions = Channel.empty()

    cram_intervals = cram.combine(intervals)
        // Move num_intervals to meta map
        .map{ meta, cram, crai, intervals, intervals_index, num_intervals -> [ meta + [ num_intervals:num_intervals ], cram, crai, intervals, intervals_index ] }

    MANTA_GERMLINE(cram_intervals, fasta, fasta_fai)

    // Figuring out if there is one or more vcf(s) from the same sample
    small_indels_vcf = MANTA_GERMLINE.out.candidate_small_indels_vcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more vcf(s) from the same sample
    sv_vcf = MANTA_GERMLINE.out.candidate_sv_vcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more vcf(s) from the same sample
    diploid_sv_vcf = MANTA_GERMLINE.out.diploid_sv_vcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    diploid_sv_vcf_to_merge = diploid_sv_vcf.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()
    small_indels_vcf_to_merge = small_indels_vcf.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()
    sv_vcf_to_merge = sv_vcf.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()

    MERGE_MANTA_DIPLOID(diploid_sv_vcf_to_merge, dict)
    MERGE_MANTA_SMALL_INDELS(small_indels_vcf_to_merge, dict)
    MERGE_MANTA_SV(sv_vcf_to_merge, dict)

    // Mix intervals and no_intervals channels together
    // Only diploid SV should get annotated
    vcf = Channel.empty().mix(MERGE_MANTA_DIPLOID.out.vcf, diploid_sv_vcf.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'manta' ], vcf ] }

    versions = versions.mix(MERGE_MANTA_DIPLOID.out.versions)
    versions = versions.mix(MERGE_MANTA_SMALL_INDELS.out.versions)
    versions = versions.mix(MERGE_MANTA_SV.out.versions)
    versions = versions.mix(MANTA_GERMLINE.out.versions)

    emit:
    vcf

    versions
}
