/*
========================================================================================
    ANNOTATION
========================================================================================
*/

params.annotation_cache             = [:]
params.bgziptabix_merge_vep_options = [:]
params.bgziptabix_snpeff_options    = [:]
params.bgziptabix_vep_options       = [:]
params.merge_vep_options            = [:]
params.snpeff_options               = [:]
params.snpeff_tag                   = [:]
params.vep_options                  = [:]
params.vep_tag                      = [:]

include { SNPEFF_ANNOTATE } from '../nf-core/snpeff_annotate' addParams(
    bgziptabix_snpeff_options: params.bgziptabix_snpeff_options,
    snpeff_options:            params.snpeff_options,
    snpeff_tag:                params.snpeff_tag,
    use_cache:                 params.annotation_cache
)

include { ENSEMBLVEP_ANNOTATE as MERGE_ANNOTATE } from '../nf-core/ensemblvep_annotate' addParams(
    bgziptabix_vep_options:    params.bgziptabix_merge_vep_options,
    use_cache:                 params.annotation_cache,
    vep_options:               params.merge_vep_options,
    vep_tag:                   params.vep_tag
)

include { ENSEMBLVEP_ANNOTATE } from '../nf-core/ensemblvep_annotate' addParams(
    bgziptabix_vep_options:    params.bgziptabix_vep_options,
    use_cache:                 params.annotation_cache,
    vep_options:               params.vep_options,
    vep_tag:                   params.vep_tag
)

workflow ANNOTATE {
    take:
    vcf          // channel: [ val(meta), vcf ]
    tools
    snpeff_db
    snpeff_cache
    vep_genome
    vep_species
    vep_cache_version
    vep_cache

    main:
    merge_vcf_ann     = Channel.empty()
    merge_vep_report  = Channel.empty()
    merge_vep_version = Channel.empty()
    snpeff_report     = Channel.empty()
    snpeff_vcf_ann    = Channel.empty()
    snpeff_version    = Channel.empty()
    vep_report        = Channel.empty()
    vep_vcf_ann       = Channel.empty()
    vep_version       = Channel.empty()

    if ('snpeff' in tools || 'merge' in tools) {
        (snpeff_vcf_ann, snpeff_report, snpeff_version) = SNPEFF_ANNOTATE(vcf, snpeff_db, snpeff_cache)
    }

    if ('merge' in tools) {
        vcf_ann_for_merge = snpeff_vcf_ann.map{ meta, vcf, tbi -> [meta, vcf] }
        (merge_vcf_ann, merge_vep_report, merge_vep_version) = MERGE_ANNOTATE(vcf_ann_for_merge, vep_genome, vep_species, vep_cache_version, vep_cache)
    }

    if ('vep' in tools) {
        (vep_vcf_ann, vep_report, vep_version) = ENSEMBLVEP_ANNOTATE(vcf, vep_genome, vep_species, vep_cache_version, vep_cache)
    }

    vcf_ann = snpeff_vcf_ann.mix(merge_vcf_ann, vep_vcf_ann)
    reports = snpeff_report.mix(merge_vep_report, vep_report)
    version = snpeff_version.first().mix(merge_vep_version.first(), vep_version.first())

    emit:
        reports
        vcf_ann
        version
}
