//
// ANNOTATION
//

include { ANNOTATION_SNPEFF                       } from '../nf-core/annotation_snpeff/main'
include { ANNOTATION_ENSEMBLVEP as MERGE_ANNOTATE } from '../nf-core/annotation_ensemblvep/main'
include { ANNOTATION_ENSEMBLVEP                   } from '../nf-core/annotation_ensemblvep/main'

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

    if (tools.contains('merge') || tools.contains('snpeff')) {
        ANNOTATION_SNPEFF(vcf, snpeff_db, snpeff_cache)

        ch_reports  = ch_reports.mix(ANNOTATION_SNPEFF.out.reports)
        ch_vcf_ann  = ch_vcf_ann.mix(ANNOTATION_SNPEFF.out.vcf_tbi)
        ch_versions = ch_versions.mix(ANNOTATION_SNPEFF.out.versions.first())
    }

    if (tools.contains('merge')) {
        vcf_ann_for_merge = ANNOTATION_SNPEFF.out.vcf_tbi.map{ meta, vcf, tbi -> [meta, vcf] }
        MERGE_ANNOTATE(vcf_ann_for_merge, vep_genome, vep_species, vep_cache_version, vep_cache)

        ch_reports  = ch_reports.mix(MERGE_ANNOTATE.out.reports)
        ch_vcf_ann  = ch_vcf_ann.mix(MERGE_ANNOTATE.out.vcf_tbi)
        ch_versions = ch_versions.mix(MERGE_ANNOTATE.out.versions.first())
    }

    if (tools.contains('vep')) {
        ANNOTATION_ENSEMBLVEP(vcf, vep_genome, vep_species, vep_cache_version, vep_cache)

        ch_reports  = ch_reports.mix(ANNOTATION_ENSEMBLVEP.out.reports)
        ch_vcf_ann  = ch_vcf_ann.mix(ANNOTATION_ENSEMBLVEP.out.vcf_tbi)
        ch_versions = ch_versions.mix(ANNOTATION_ENSEMBLVEP.out.versions.first())
    }

    emit:
    vcf_ann  = ch_vcf_ann      // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    reports  = ch_reports      //    path: *.html
    versions = ch_versions     //    path: versions.yml
}
