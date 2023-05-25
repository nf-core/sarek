include { GATK4_MERGEVCFS as MERGE_DEEPVARIANT_GVCF } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_DEEPVARIANT_VCF  } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { DEEPVARIANT                               } from '../../../modules/nf-core/deepvariant/main'
include { TABIX_TABIX as TABIX_VC_DEEPVARIANT_GVCF  } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_VC_DEEPVARIANT_VCF   } from '../../../modules/nf-core/tabix/tabix/main'

// Deepvariant: https://github.com/google/deepvariant/issues/510
workflow BAM_VARIANT_CALLING_DEEPVARIANT {
    take:
    cram          // channel: [mandatory] [ meta, cram, crai ]
    dict          // channel: [optional]  [ meta, dict ]
    fasta         // channel: [mandatory] [ fasta ]
    fasta_fai     // channel: [mandatory] [ fasta_fai ]
    intervals     // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals

    main:
    versions = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    cram_intervals = cram.combine(intervals)
        // Move num_intervals to meta map
        .map{ meta, cram, crai, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], cram, crai, intervals ]}

    DEEPVARIANT(cram_intervals, fasta.map{ fasta -> [ [ id:fasta.baseName ], fasta ] }, fasta_fai.map{ fasta_fai -> [ [ id:fasta_fai.baseName ], fasta_fai ] }, [ [ id:'null' ], [] ])

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_out = DEEPVARIANT.out.vcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more gvcf(s) from the same sample
    gvcf_out = DEEPVARIANT.out.gvcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    gvcf_to_merge = gvcf_out.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()
    vcf_to_merge = vcf_out.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()

    MERGE_DEEPVARIANT_GVCF(gvcf_to_merge, dict)
    MERGE_DEEPVARIANT_VCF(vcf_to_merge, dict)

    // Only when no_intervals
    TABIX_VC_DEEPVARIANT_GVCF(gvcf_out.no_intervals)
    TABIX_VC_DEEPVARIANT_VCF(vcf_out.no_intervals)

    // Mix intervals and no_intervals channels together
    gvcf = Channel.empty().mix(MERGE_DEEPVARIANT_GVCF.out.vcf, gvcf_out.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'deepvariant' ], vcf ] }

    // Mix intervals and no_intervals channels together
    vcf = Channel.empty().mix(MERGE_DEEPVARIANT_VCF.out.vcf, vcf_out.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'deepvariant' ], vcf ] }

    versions = versions.mix(DEEPVARIANT.out.versions)
    versions = versions.mix(MERGE_DEEPVARIANT_GVCF.out.versions)
    versions = versions.mix(MERGE_DEEPVARIANT_VCF.out.versions)
    versions = versions.mix(TABIX_VC_DEEPVARIANT_GVCF.out.versions)
    versions = versions.mix(TABIX_VC_DEEPVARIANT_VCF.out.versions)

    emit:
    gvcf
    vcf

    versions
}
