include { CAT_CAT as CAT_MPILEUP         } from '../../../modules/nf-core/cat/cat/main'
include { BCFTOOLS_MPILEUP               } from '../../../modules/nf-core/bcftools/mpileup/main'
include { SAMTOOLS_MPILEUP               } from '../../../modules/nf-core/samtools/mpileup/main'
include { GATK4_MERGEVCFS                } from '../../../modules/nf-core/gatk4/mergevcfs/main'


workflow BAM_VARIANT_CALLING_MPILEUP {
    take:
    cram      // channel: [mandatory] [ meta, cram, crai ]
    dict      // channel: [mandatory] [ meta, dict ]
    fasta     // channel: [mandatory] [ fasta ]
    intervals // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals

    main:
    versions = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    cram_intervals = cram.combine(intervals)
        // Move num_intervals to meta map and reorganize channel for BCFTOOLS_MPILEUP/SAMTOOLS_MPILEUP modules
        .map{ meta, cram, crai, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], cram, intervals ] }

    keep_bcftools_mpileup = false
    BCFTOOLS_MPILEUP(cram_intervals, fasta, keep_bcftools_mpileup)

    SAMTOOLS_MPILEUP(cram_intervals, fasta)

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_mpileup = BCFTOOLS_MPILEUP.out.vcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more mpileup(s) from the same sample
    mpileup_samtools = SAMTOOLS_MPILEUP.out.mpileup.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Merge mpileup and natural order sort them
    mpileup_to_merge = mpileup_samtools.intervals.map{ meta, pileup -> [ groupKey(meta, meta.num_intervals), pileup ] }.groupTuple(sort:true)
    CAT_MPILEUP(mpileup_to_merge)

    // Merge VCF
    vcf_to_merge = vcf_mpileup.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ] }.groupTuple()
    GATK4_MERGEVCFS(vcf_to_merge, dict)

    // Mix intervals and no_intervals channels together
    mpileup = CAT_MPILEUP.out.file_out.mix(mpileup_samtools.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, mpileup -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'mpileup' ], mpileup ] }
    vcf = GATK4_MERGEVCFS.out.vcf.mix(vcf_mpileup.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'mpileup' ], vcf ] }

    versions = versions.mix(SAMTOOLS_MPILEUP.out.versions)
    versions = versions.mix(BCFTOOLS_MPILEUP.out.versions)
    versions = versions.mix(CAT_MPILEUP.out.versions)
    versions = versions.mix(GATK4_MERGEVCFS.out.versions)

    emit:
    mpileup
    vcf

    versions
}
