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

include { SNPEFF_ANNOTATE } from '../nf-core/snpeff' addParams(
    bgziptabix_snpeff_options: params.bgziptabix_snpeff_options,
    snpeff_options:            params.snpeff_options,
    snpeff_tag:                params.snpeff_tag,
    use_cache:                 params.annotation_cache
)

include { VEP_ANNOTATE as MERGE_ANNOTATE } from '../nf-core/vep' addParams(
    bgziptabix_vep_options:    params.bgziptabix_merge_vep_options,
    use_cache:                 params.annotation_cache,
    vep_options:               params.merge_vep_options,
    vep_tag:                   params.vep_tag
)

include { VEP_ANNOTATE } from '../nf-core/vep' addParams(
    bgziptabix_vep_options:    params.bgziptabix_vep_options,
    use_cache:                 params.annotation_cache,
    vep_options:               params.vep_options,
    vep_tag:                   params.vep_tag
)

workflow ANNOTATE {
    take:
    vcf          // channel: [ val(meta), vcf, tbi ]
    tools
    snpeff_db
    snpeff_cache
    vep_genome
    vep_species
    vep_cache_version
    vep_cache

    main:
    snpeff_vcf_ann = Channel.empty()
    merge_vcf_ann = Channel.empty()
    vep_vcf_ann = Channel.empty()

    if ('snpeff' in tools || 'merge' in tools) {
        (snpeff_vcf_ann, snpeff_reports, snpeff_version) = SNPEFF_ANNOTATE(vcf, snpeff_db, snpeff_cache)
    }

    if ('merge' in tools) {
        vcf_ann_for_merge = snpeff_vcf_ann.map{ meta, vcf, tbi -> [meta, vcf] } 
        (merge_vcf_ann, merge_vep_reports, merge_vep_version) = MERGE_ANNOTATE(vcf_ann_for_merge, vep_genome, vep_species, vep_cache_version, vep_cache)
    }

    if ('vep' in tools) {
        (vep_vcf_ann, vep_reports, vep_version) = VEP_ANNOTATE(vcf, vep_genome, vep_species, vep_cache_version, vep_cache)
    }

    vcf_ann = snpeff_vcf_ann.mix(merge_vcf_ann, vep_vcf_ann)

    emit:
        vcf_ann
}
