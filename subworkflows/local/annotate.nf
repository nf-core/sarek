//
// ANNOTATION
//

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
    ch_reports  = Channel.empty()
    ch_vcf_ann  = Channel.empty()
    ch_versions = Channel.empty()

    if ('snpeff' in tools || 'merge' in tools) {
        SNPEFF_ANNOTATE(vcf, snpeff_db, snpeff_cache)

        ch_reports  = ch_reports.mix(SNPEFF_ANNOTATE.out.reports)
        ch_vcf_ann  = ch_vcf_ann.mix(SNPEFF_ANNOTATE.out.vcf_tbi)
        ch_versions = ch_versions.mix(SNPEFF_ANNOTATE.out.versions.first())
    }

    if ('merge' in tools) {
        vcf_ann_for_merge = SNPEFF_ANNOTATE.out.vcf_tbi.map{ meta, vcf, tbi -> [meta, vcf] }
        MERGE_ANNOTATE(vcf_ann_for_merge, vep_genome, vep_species, vep_cache_version, vep_cache)

        ch_reports  = ch_reports.mix(MERGE_ANNOTATE.out.reports)
        ch_vcf_ann  = ch_vcf_ann.mix(MERGE_ANNOTATE.out.vcf_tbi)
        ch_versions = ch_versions.mix(MERGE_ANNOTATE.out.versions.first())
    }

    if ('vep' in tools) {
        ENSEMBLVEP_ANNOTATE(vcf, vep_genome, vep_species, vep_cache_version, vep_cache)

        ch_reports  = ch_reports.mix(ENSEMBLVEP_ANNOTATE.out.reports)
        ch_vcf_ann  = ch_vcf_ann.mix(ENSEMBLVEP_ANNOTATE.out.vcf_tbi)
        ch_versions = ch_versions.mix(ENSEMBLVEP_ANNOTATE.out.versions.first())
    }

    emit:
    vcf_ann  = ch_vcf_ann      // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    reports  = ch_reports      //    path: *.html
    versions = ch_versions     //    path: versions.yml
}
