include { GATK4_MERGEVCFS as MERGE_STRELKA        } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_STRELKA_GENOME } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { STRELKA_GERMLINE as STRELKA_SINGLE      } from '../../../modules/nf-core/strelka/germline/main'

workflow BAM_VARIANT_CALLING_SINGLE_STRELKA {
    take:
    cram          // channel: [mandatory] [ meta, cram, crai ]
    dict          // channel: [optional]  [ meta, dict ]
    fasta         // channel: [mandatory] [ fasta ]
    fasta_fai     // channel: [mandatory] [ fasta_fai ]
    intervals     // channel: [mandatory] [ interval.bed.gz, interval.bed.gz.tbi, num_intervals ] or [ [], [], 0 ] if no intervals

    main:
    versions = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    cram_intervals = cram.combine(intervals)
        // Move num_intervals to meta map
        .map{ meta, cram, crai, intervals, intervals_index, num_intervals -> [ meta + [ num_intervals:num_intervals ], cram, crai, intervals, intervals_index ] }

    STRELKA_SINGLE(cram_intervals, fasta, fasta_fai)

    // Figuring out if there is one or more vcf(s) from the same sample
    genome_vcf = STRELKA_SINGLE.out.genome_vcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf = STRELKA_SINGLE.out.vcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    genome_vcf_to_merge = genome_vcf.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()
    vcf_to_merge = vcf.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()

    MERGE_STRELKA(vcf_to_merge, dict)
    MERGE_STRELKA_GENOME(genome_vcf_to_merge, dict)

    // Mix intervals and no_intervals channels together
    // Only strelka variant vcf SV should get annotated
    vcf = Channel.empty().mix(MERGE_STRELKA.out.vcf, vcf.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'strelka' ], vcf ] }

    versions = versions.mix(MERGE_STRELKA.out.versions)
    versions = versions.mix(MERGE_STRELKA_GENOME.out.versions)
    versions = versions.mix(STRELKA_SINGLE.out.versions)

    emit:
    vcf

    versions
}
