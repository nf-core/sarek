//
// FREEBAYES variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BCFTOOLS_SORT                         } from '../../../modules/nf-core/bcftools/sort/main'
include { FREEBAYES                             } from '../../../modules/nf-core/freebayes/main'
include { GATK4_MERGEVCFS as MERGE_FREEBAYES    } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { TABIX_TABIX     as TABIX_VC_FREEBAYES } from '../../../modules/nf-core/tabix/tabix/main'

workflow BAM_VARIANT_CALLING_FREEBAYES {
    take:
    cram      // channel: [mandatory] [ meta, cram1, crai1, cram2, crai2 ] or [ meta, cram, crai, [], [] ]
    dict      // channel: [mandatory] [ meta, dict ]
    fasta     // channel: [mandatory] [ fasta ]
    fasta_fai // channel: [mandatory] [ fasta_fai ]
    intervals // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals

    main:
    versions = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    cram_intervals = cram.combine(intervals)
        // Move num_intervals to meta map and reorganize channel for FREEBAYES module
        .map{ meta, cram1, crai1, cram2, crai2, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], cram1, crai1, cram2, crai2, intervals ]}

    FREEBAYES(cram_intervals, fasta, fasta_fai, [], [], [])

    BCFTOOLS_SORT(FREEBAYES.out.vcf)

    // Figuring out if there is one or more vcf(s) from the same sample
    bcftools_vcf_out = BCFTOOLS_SORT.out.vcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    vcf_to_merge = bcftools_vcf_out.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()
    MERGE_FREEBAYES(vcf_to_merge, dict)

    // Only when no_intervals
    TABIX_VC_FREEBAYES(bcftools_vcf_out.no_intervals)

    // Mix intervals and no_intervals channels together
    vcf = MERGE_FREEBAYES.out.vcf.mix(bcftools_vcf_out.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'freebayes' ], vcf ] }

    versions = versions.mix(BCFTOOLS_SORT.out.versions)
    versions = versions.mix(MERGE_FREEBAYES.out.versions)
    versions = versions.mix(FREEBAYES.out.versions)
    versions = versions.mix(TABIX_VC_FREEBAYES.out.versions)

    emit:
    vcf

    versions
}
