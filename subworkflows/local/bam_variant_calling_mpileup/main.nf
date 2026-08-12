//
// MPILEUP variant calling: BCFTOOLS for variantcalling, SAMTools for controlfreec input
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BCFTOOLS_MPILEUP                          } from '../../../modules/nf-core/bcftools/mpileup'
include { FIND_CONCATENATE  as CAT_MPILEUP          } from '../../../modules/nf-core/find/concatenate'
include { GATK4_MERGEVCFS as MERGE_BCFTOOLS_MPILEUP } from '../../../modules/nf-core/gatk4/mergevcfs'
include { SAMTOOLS_MPILEUP                          } from '../../../modules/nf-core/samtools/mpileup'

workflow BAM_VARIANT_CALLING_MPILEUP {
    take:
    cram      // channel: [mandatory] [ meta, cram, crai ]
    dict      // channel: [mandatory] [ meta, dict ]
    fasta     // channel: [mandatory] [ meta, fasta ]
    fai       // channel: [mandatory] [ meta, fasta_fai ]
    intervals // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals

    main:

    // Combine cram and intervals for spread and gather strategy
    cram_intervals = cram
        .combine(intervals)
        .map { meta, cram_, _crai, intervals_, num_intervals -> [meta + [num_intervals: num_intervals], cram_, intervals_] }

    // per-chunk bed restricts mpileup only; intervals_call left empty
    cram_intervals_bcftools = cram_intervals.map { meta, cram_, intervals_ -> [meta, cram_, intervals_, []] }

    // samtools/mpileup now expects [ meta, input, index, intervals ]; keep the crai
    cram_crai_intervals = cram.combine(intervals).map { meta, cram_, crai, intervals_, num_intervals -> [meta + [num_intervals: num_intervals], cram_, crai, intervals_] }

    fasta_fai = fasta.combine(fai).map { meta, fasta_, _meta2, fai_ -> [meta, fasta_, fai_] }.collect()

    // Run, if --tools mpileup
    keep_bcftools_mpileup = false
    BCFTOOLS_MPILEUP(cram_intervals_bcftools, fasta_fai, keep_bcftools_mpileup)

    //Only run, if --tools ControlFreec
    SAMTOOLS_MPILEUP(cram_crai_intervals, fasta_fai)

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_mpileup = BCFTOOLS_MPILEUP.out.vcf.branch { meta, _vcf ->
        intervals: meta.num_intervals > 1
        no_intervals: meta.num_intervals <= 1
    }

    // Figuring out if there is one or more tbi(s) from the same sample
    tbi_mpileup = BCFTOOLS_MPILEUP.out.index.branch { meta, _tbi ->
        intervals: meta.num_intervals > 1
        no_intervals: meta.num_intervals <= 1
    }

    // Figuring out if there is one or more mpileup(s) from the same sample
    mpileup_samtools = SAMTOOLS_MPILEUP.out.mpileup.branch { meta, _mpileup ->
        intervals: meta.num_intervals > 1
        no_intervals: meta.num_intervals <= 1
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
    tbi = MERGE_BCFTOOLS_MPILEUP.out.tbi
        .mix(tbi_mpileup.no_intervals)
        .map { meta, tbi -> [meta - meta.subMap('num_intervals') + [variantcaller: 'bcftools'], tbi] }

    emit:
    mpileup
    vcf
    tbi
}
