//
// JOINT GERMLINE CALLING for DeepVariant
//
// Merge per-sample DeepVariant gVCFs and perform joint genotyping with GLnexus,
// then convert the cohort BCF to a bgzipped, indexed multi-sample VCF.
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GLNEXUS                                } from '../../../modules/nf-core/glnexus/main'
include { BCFTOOLS_VIEW as GLNEXUS_BCFTOOLS_VIEW } from '../../../modules/nf-core/bcftools/view/main'

workflow BAM_JOINT_GERMLINE_DEEPVARIANT {
    take:
    gvcf                   // channel: [mandatory] [ meta, gvcf ] one entry per sample
    intervals_bed_combined // channel: [mandatory] [ intervals.bed ] or [] if no_intervals

    main:
    versions = Channel.empty()

    // Gather all per-sample gVCFs into a single cohort call for GLnexus
    // Sort by sample id so the GLnexus input order is deterministic (reproducible output)
    // The trailing [] is the optional GLnexus custom_config (the model is selected via ext.args)
    joint_gvcfs = gvcf
        .map { meta, vcf -> [ meta.sample, vcf ] }
        .toSortedList { a, b -> a[0] <=> b[0] }
        .map { samples ->
            if (samples.size() < 2) {
                error("Joint genotyping requires at least 2 DeepVariant gVCFs, but got ${samples.size()}")
            }
            [ [ id:'joint_variant_calling' ], samples.collect { it[1] }, [] ]
        }

    // Restrict GLnexus to the combined calling intervals (single combined bed -> path, matching path(bed))
    bed = intervals_bed_combined.map { it -> it ? [ [ id:it[0].baseName ], it[0] ] : [ [ id:'no_intervals' ], [] ] }

    // Joint genotyping: per-sample gVCFs -> one multi-sample cohort BCF
    GLNEXUS(joint_gvcfs, bed)

    // Convert the cohort BCF to a bgzipped, indexed multi-sample VCF
    // GLnexus output is already position-sorted, so a plain view (+ index) is enough
    // bcftools view takes [ meta, vcf, index ]; the BCF has no index, so pass []
    GLNEXUS_BCFTOOLS_VIEW(
        GLNEXUS.out.bcf.map { meta, bcf -> [ meta, bcf, [] ] },
        [],
        [],
        [],
    )

    // Rework meta for variantscalled.csv and downstream annotation tools
    genotype_vcf = GLNEXUS_BCFTOOLS_VIEW.out.vcf.map { _meta, vcf ->
        [ [ id:'joint_variant_calling', patient:'all_samples', variantcaller:'deepvariant' ], vcf ]
    }
    genotype_index = GLNEXUS_BCFTOOLS_VIEW.out.tbi.map { _meta, tbi ->
        [ [ id:'joint_variant_calling', patient:'all_samples', variantcaller:'deepvariant' ], tbi ]
    }

    // GLnexus emits its version through the global "versions" topic channel
    versions = versions.mix(GLNEXUS_BCFTOOLS_VIEW.out.versions)

    emit:
    genotype_index // channel: [ val(meta), tbi ]
    genotype_vcf   // channel: [ val(meta), vcf ]

    versions       // channel: [ versions.yml ]
}
