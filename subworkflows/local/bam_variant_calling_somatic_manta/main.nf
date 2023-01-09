include { GATK4_MERGEVCFS as MERGE_MANTA_DIPLOID           } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_MANTA_SMALL_INDELS      } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_MANTA_SOMATIC           } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_MANTA_SV                } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { MANTA_SOMATIC                                    } from '../../../modules/nf-core/manta/somatic/main'

workflow BAM_VARIANT_CALLING_SOMATIC_MANTA {
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
        .map{ meta, cram, crai, intervals, intervals_index, num_intervals -> [ meta + [ num_intervals:num_intervals ], cram, crai, intervals, intervals_index ]}

    MANTA_SOMATIC(cram_intervals, fasta, fasta_fai)

    // Figuring out if there is one or more vcf(s) from the same sample
    candidate_small_indels_vcf = MANTA_SOMATIC.out.candidate_small_indels_vcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more vcf(s) from the same sample
    candidate_small_indels_vcf_tbi = MANTA_SOMATIC.out.candidate_small_indels_vcf_tbi.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more vcf(s) from the same sample
    candidate_sv_vcf = MANTA_SOMATIC.out.candidate_sv_vcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more vcf(s) from the same sample
    diploid_sv_vcf = MANTA_SOMATIC.out.diploid_sv_vcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more vcf(s) from the same sample
    somatic_sv_vcf = MANTA_SOMATIC.out.somatic_sv_vcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    candidate_small_indels_vcf_to_merge = candidate_small_indels_vcf.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()
    candidate_sv_vcf_to_merge = candidate_sv_vcf.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()
    diploid_sv_vcf_to_merge = diploid_sv_vcf.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()
    somatic_sv_vcf_to_merge = somatic_sv_vcf.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()

    MERGE_MANTA_SMALL_INDELS(candidate_small_indels_vcf_to_merge, dict)
    MERGE_MANTA_SV(candidate_sv_vcf_to_merge, dict)
    MERGE_MANTA_DIPLOID(diploid_sv_vcf_to_merge, dict)
    MERGE_MANTA_SOMATIC(somatic_sv_vcf_to_merge, dict)

    // Mix intervals and no_intervals channels together
    // Only diploid and somatic SV should get annotated
    vcf = Channel.empty().mix(MERGE_MANTA_DIPLOID.out.vcf, MERGE_MANTA_SOMATIC.out.vcf, diploid_sv_vcf.no_intervals, somatic_sv_vcf.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'manta' ], vcf ] }

    // Mix intervals and no_intervals channels together
    // Only joining reads for StrelkaBP
    candidate_small_indels_vcf = Channel.empty().mix(MERGE_MANTA_SMALL_INDELS.out.vcf, candidate_small_indels_vcf.no_intervals)
        // remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals'), vcf ] }

    // Mix intervals and no_intervals channels together
    // Only joining reads for StrelkaBP
    candidate_small_indels_vcf_tbi = Channel.empty().mix(MERGE_MANTA_SMALL_INDELS.out.tbi, candidate_small_indels_vcf_tbi.no_intervals)
        // remove no longer necessary field: num_intervals
        .map{ meta, vcf, tbi -> [ meta - meta.subMap('num_intervals'), vcf, tbi ] }

    versions = versions.mix(MERGE_MANTA_SV.out.versions)
    versions = versions.mix(MERGE_MANTA_SMALL_INDELS.out.versions)
    versions = versions.mix(MERGE_MANTA_DIPLOID.out.versions)
    versions = versions.mix(MERGE_MANTA_SOMATIC.out.versions)
    versions = versions.mix(MANTA_SOMATIC.out.versions)

    emit:
    candidate_small_indels_vcf
    candidate_small_indels_vcf_tbi
    vcf

    versions
}
