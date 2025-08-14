//
// FREEBAYES variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BCFTOOLS_SORT                         } from '../../../modules/nf-core/bcftools/sort'
include { FREEBAYES                             } from '../../../modules/nf-core/freebayes'
include { GATK4_MERGEVCFS as MERGE_FREEBAYES    } from '../../../modules/nf-core/gatk4/mergevcfs'
include { TABIX_TABIX     as TABIX_VC_FREEBAYES } from '../../../modules/nf-core/tabix/tabix'
include { VCFLIB_VCFFILTER                      } from '../../../modules/nf-core/vcflib/vcffilter'

workflow BAM_VARIANT_CALLING_FREEBAYES {
    take:
    ch_cram      // channel: [mandatory] [ meta, cram1, crai1, cram2, crai2 ] or [ meta, cram, crai, [], [] ]
    ch_dict      // channel: [mandatory] [ meta, dict ]
    ch_fasta     // channel: [mandatory] [ meta, fasta ]
    ch_fasta_fai // channel: [mandatory] [ meta, fasta_fai ]
    ch_intervals // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals

    main:
    ch_versions= Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    cram_intervals = ch_cram.combine(ch_intervals)
        // Move num_intervals to meta map and reorganize channel for FREEBAYES module
        .map{ meta, cram1, crai1, cram2, crai2, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], cram1, crai1, cram2, crai2, intervals ]}

    FREEBAYES(cram_intervals, ch_fasta, ch_fasta_fai, [[id:'null'], []], [[id:'null'], []], [[id:'null'], []])

    BCFTOOLS_SORT(FREEBAYES.out.vcf)

    // Figuring out if there is one or more vcf(s) from the same sample
    bcftools_vcf_out = BCFTOOLS_SORT.out.vcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    vcf_to_merge = bcftools_vcf_out.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()
    MERGE_FREEBAYES(vcf_to_merge, ch_dict)

    // Only when no_intervals
    TABIX_VC_FREEBAYES(bcftools_vcf_out.no_intervals)

    // Mix intervals and no_intervals channels together
    ch_vcf = MERGE_FREEBAYES.out.vcf.mix(bcftools_vcf_out.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'freebayes' ], vcf ] }

    vcf_filtered = VCFLIB_VCFFILTER(ch_vcf.map{ meta, vcf -> [ meta, vcf, [] ] })

    ch_versions= ch_versions.mix(BCFTOOLS_SORT.out.versions)
    ch_versions= ch_versions.mix(MERGE_FREEBAYES.out.versions)
    ch_versions= ch_versions.mix(FREEBAYES.out.versions)
    ch_versions= ch_versions.mix(TABIX_VC_FREEBAYES.out.versions)
    ch_versions= ch_versions.mix(VCFLIB_VCFFILTER.out.versions)

    emit:
    ch_vcf
    vcf_filtered // channel: [ meta, vcf ]
    ch_versions
}
