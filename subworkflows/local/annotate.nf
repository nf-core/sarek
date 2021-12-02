//
// ANNOTATION
//

include { SNPEFF_ANNOTATE                       } from '../nf-core/snpeff_annotate'
include { ENSEMBLVEP_ANNOTATE as MERGE_ANNOTATE } from '../nf-core/ensemblvep_annotate'
include { ENSEMBLVEP_ANNOTATE                   } from '../nf-core/ensemblvep_annotate'

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
        SNPEFF_ANNOTATE(vcf, snpeff_db, snpeff_cache)

        ch_reports  = ch_reports.mix(SNPEFF_ANNOTATE.out.reports)
        ch_vcf_ann  = ch_vcf_ann.mix(SNPEFF_ANNOTATE.out.vcf_tbi)
        ch_versions = ch_versions.mix(SNPEFF_ANNOTATE.out.versions.first())
    }

    if (tools.contains('merge')) {
        vcf_ann_for_merge = SNPEFF_ANNOTATE.out.vcf_tbi.map{ meta, vcf, tbi -> [meta, vcf] }
        MERGE_ANNOTATE(vcf_ann_for_merge, vep_genome, vep_species, vep_cache_version, vep_cache)

        ch_reports  = ch_reports.mix(MERGE_ANNOTATE.out.reports)
        ch_vcf_ann  = ch_vcf_ann.mix(MERGE_ANNOTATE.out.vcf_tbi)
        ch_versions = ch_versions.mix(MERGE_ANNOTATE.out.versions.first())
    }

    if (tools.contains('vep')) {
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
