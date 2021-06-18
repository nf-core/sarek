/*
========================================================================================
    ANNOTATION
========================================================================================
*/

params.snpeff_options               = [:]
params.vep_options                  = [:]
params.vep_merge_options            = [:]
params.bgziptabix_snpeff_options    = [:]
params.bgziptabix_vep_options       = [:]
params.bgziptabix_vep_merge_options = [:]

include { SNPEFF_ANNOTATE } from '../nf-core/snpeff' addParams(
    snpeff_options:            params.snpeff_options,
    bgziptabix_snpeff_options: params.bgziptabix_snpeff_options
)

include { VEP_ANNOTATE }    from '../nf-core/vep' addParams(
    vep_options:               params.vep_options,
    bgziptabix_vep_options:    params.bgziptabix_vep_options
)

include { VEP_ANNOTATE as VEP_MERGE_ANNOTATE }    from '../nf-core/vep' addParams(
    vep_merge_options:               params.vep_merge_options,
    bgziptabix_vep_merge_options:    params.bgziptabix_vep_merge_options
)

workflow ANNOTATE {
    take:
    vcf          // channel: [ val(meta), vcf, tbi ]
    use_cache    //    bool: use cache (default: false)
    tools
    snpeff_db
    snpeff_cache
    snpeff_tag
    vep_genome
    vep_species
    vep_cache_version
    vep_cache
    vep_tag

    main:
    snpeff_vcf_ann = Channel.empty()
    merge_vcf_ann = Channel.empty()
    vep_vcf_ann = Channel.empty()

    if ('snpeff' in tools || 'merge' in tools) {
        (snpeff_vcf_ann, snpeff_reports, snpeff_version) = SNPEFF_ANNOTATE(vcf, snpeff_db, use_cache, snpeff_cache, snpeff_tag)
    }

    if ('merge' in tools) {
        vcf_ann_for_merge = snpeff_vcf_ann.map{ meta, vcf, tbi -> [meta, vcf] } 
        (merge_vcf_ann, vep_reports, vep_version) = VEP_MERGE_ANNOTATE(vcf_ann_for_merge, vep_genome, vep_species, vep_cache_version, use_cache, vep_cache, vep_tag)
    }

    if ('vep' in tools) {
        (vep_vcf_ann, vep_reports, vep_version) = VEP_ANNOTATE(vcf, vep_genome, vep_species, vep_cache_version, use_cache, vep_cache, vep_tag)
    }

    vcf_ann = snpeff_vcf_ann.mix(merge_vcf_ann, vep_vcf_ann)

    emit:
        vcf_ann
}
