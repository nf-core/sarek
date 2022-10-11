//
// ANNOTATION
//

include { ANNOTATION_SNPEFF                         } from '../nf-core/annotation/snpeff/main'
include { ANNOTATION_ENSEMBLVEP as ANNOTATION_MERGE } from '../nf-core/annotation/ensemblvep/main'
include { ANNOTATION_ENSEMBLVEP                     } from '../nf-core/annotation/ensemblvep/main'

workflow VCF_ANNOTATE_ALL {
    take:
    vcf          // channel: [ val(meta), vcf ]
    fasta
    tools        // Mandatory, list of tools to apply
    snpeff_db
    snpeff_cache
    vep_genome
    vep_species
    vep_cache_version
    vep_cache
    vep_extra_files

    main:
    ch_reports   = Channel.empty()
    ch_vcf_ann   = Channel.empty()
    ch_tab_ann   = Channel.empty()
    ch_json_ann  = Channel.empty()
    ch_versions  = Channel.empty()

    if (tools.split(',').contains('merge') || tools.split(',').contains('snpeff')) {
        ANNOTATION_SNPEFF(vcf, snpeff_db, snpeff_cache)

        ch_reports  = ch_reports.mix(ANNOTATION_SNPEFF.out.reports)
        ch_vcf_ann  = ch_vcf_ann.mix(ANNOTATION_SNPEFF.out.vcf_tbi)
        ch_versions = ch_versions.mix(ANNOTATION_SNPEFF.out.versions.first())
    }

    if (tools.split(',').contains('merge')) {
        vcf_ann_for_merge = ANNOTATION_SNPEFF.out.vcf_tbi.map{ meta, vcf, tbi -> [meta, vcf] }
        ANNOTATION_MERGE(vcf_ann_for_merge, fasta, vep_genome, vep_species, vep_cache_version, vep_cache, vep_extra_files)

        ch_reports  = ch_reports.mix(ANNOTATION_MERGE.out.reports)
        ch_vcf_ann  = ch_vcf_ann.mix(ANNOTATION_MERGE.out.vcf_tbi)
        ch_versions = ch_versions.mix(ANNOTATION_MERGE.out.versions.first())
    }

    if (tools.split(',').contains('vep')) {
        ANNOTATION_ENSEMBLVEP(vcf, fasta, vep_genome, vep_species, vep_cache_version, vep_cache, vep_extra_files)

        ch_reports   = ch_reports.mix(ANNOTATION_ENSEMBLVEP.out.reports)
        ch_vcf_ann   = ch_vcf_ann.mix(ANNOTATION_ENSEMBLVEP.out.vcf_tbi)
        ch_tab_ann   = ch_vcf_ann.mix(ANNOTATION_ENSEMBLVEP.out.tab)
        ch_json_ann  = ch_vcf_ann.mix(ANNOTATION_ENSEMBLVEP.out.json)
        ch_versions  = ch_versions.mix(ANNOTATION_ENSEMBLVEP.out.versions.first())
    }

    emit:
    vcf_ann   = ch_vcf_ann      // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    tab_ann   = ch_tab_ann
    json_ann  = ch_json_ann
    reports   = ch_reports      //    path: *.html
    versions  = ch_versions     //    path: versions.yml
}
