//
// MPILEUP variant calling: BCFTOOLS for variantcalling, SAMTools for controlfreec input
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BCFTOOLS_MPILEUP                          } from '../../../modules/nf-core/bcftools/mpileup'
include { CAT_CAT as CAT_MPILEUP                    } from '../../../modules/nf-core/cat/cat'
include { GATK4_MERGEVCFS as MERGE_BCFTOOLS_MPILEUP } from '../../../modules/nf-core/gatk4/mergevcfs'
include { SAMTOOLS_MPILEUP                          } from '../../../modules/nf-core/samtools/mpileup'

workflow BAM_VARIANT_CALLING_MPILEUP {
    take:
    cram      // channel: [mandatory] [ meta, cram, crai ]
    dict      // channel: [mandatory] [ meta, dict ]
    fasta     // channel: [mandatory] [ fasta ]
    intervals // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals

    main:
    versions = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    cram_intervals = cram
        .combine(intervals)
        .map { meta, cram_, _crai, intervals_, num_intervals -> [meta + [num_intervals: num_intervals], cram_, intervals_] }

    // Run, if --tools mpileup
    keep_bcftools_mpileup = false
    BCFTOOLS_MPILEUP(cram_intervals, fasta, keep_bcftools_mpileup)

    //Only run, if --tools ControlFreec
    SAMTOOLS_MPILEUP(cram_intervals, fasta)

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_mpileup = BCFTOOLS_MPILEUP.out.vcf.branch {
        intervals: it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more mpileup(s) from the same sample
    mpileup_samtools = SAMTOOLS_MPILEUP.out.mpileup.branch {
        intervals: it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Merge mpileup and natural order sort them
    mpileup_to_merge = mpileup_samtools.intervals.map { meta, pileup -> [groupKey(meta, meta.num_intervals), pileup] }.groupTuple(sort: true)
    CAT_MPILEUP(mpileup_to_merge)

    // Merge VCF
    vcf_to_merge = vcf_mpileup.intervals.map { meta, vcf -> [groupKey(meta, meta.num_intervals), vcf] }.groupTuple()
    MERGE_BCFTOOLS_MPILEUP(vcf_to_merge, dict)

    // Mix intervals and no_intervals channels together
    mpileup = CAT_MPILEUP.out.file_out
        .mix(mpileup_samtools.no_intervals)
        .map { meta, mpileup -> [meta - meta.subMap('num_intervals') + [variantcaller: 'samtools'], mpileup] }
    vcf = MERGE_BCFTOOLS_MPILEUP.out.vcf
        .mix(vcf_mpileup.no_intervals)
        .map { meta, vcf -> [meta - meta.subMap('num_intervals') + [variantcaller: 'bcftools'], vcf] }

    versions = versions.mix(SAMTOOLS_MPILEUP.out.versions)
    versions = versions.mix(BCFTOOLS_MPILEUP.out.versions)
    versions = versions.mix(CAT_MPILEUP.out.versions)
    versions = versions.mix(MERGE_BCFTOOLS_MPILEUP.out.versions)

    emit:
    mpileup
    vcf
    versions
}
