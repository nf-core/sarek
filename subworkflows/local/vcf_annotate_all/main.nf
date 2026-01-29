//
// ANNOTATION
//

include { BCFTOOLS_ANNOTATE                             } from '../../../modules/nf-core/bcftools/annotate'
include { VCF_ANNOTATE_ENSEMBLVEP                       } from '../../nf-core/vcf_annotate_ensemblvep'
include { VCF_ANNOTATE_ENSEMBLVEP as VCF_ANNOTATE_MERGE } from '../../nf-core/vcf_annotate_ensemblvep'
include { VCF_ANNOTATE_SNPEFF                           } from '../../nf-core/vcf_annotate_snpeff'
include { VCF_ANNOTATE_SNPSIFT                          } from '../vcf_annotate_snpsift/main'

workflow VCF_ANNOTATE_ALL {
    take:
    vcf                        // channel: [ val(meta), vcf ]
    fasta
    tools                      // Mandatory, list of tools to apply
    snpeff_db
    snpeff_cache
    vep_genome
    vep_species
    vep_cache_version
    vep_cache
    vep_extra_files
    bcftools_annotations
    bcftools_annotations_index
    bcftools_columns
    bcftools_header_lines
    snpsift_db_configs

    main:
    reports = Channel.empty()
    vcf_ann = Channel.empty()
    tab_ann = Channel.empty()
    json_ann = Channel.empty()
    versions = Channel.empty()

    if (tools.split(',').contains('bcfann')) {
        BCFTOOLS_ANNOTATE(
            vcf.map { meta, vcf_ -> [meta, vcf_, []] }.combine(bcftools_annotations).combine(bcftools_annotations_index),
            bcftools_columns,
            bcftools_header_lines,
            [],
        )

        vcf_ann = vcf_ann.mix(BCFTOOLS_ANNOTATE.out.vcf.join(BCFTOOLS_ANNOTATE.out.tbi, failOnDuplicate: true, failOnMismatch: true))
        versions = versions.mix(BCFTOOLS_ANNOTATE.out.versions)
    }

    if (tools.split(',').contains('merge') || tools.split(',').contains('snpeff')) {
        VCF_ANNOTATE_SNPEFF(vcf, snpeff_db, snpeff_cache)

        reports = reports.mix(VCF_ANNOTATE_SNPEFF.out.reports.map { _meta, reports_ -> [reports_] })
        vcf_ann = vcf_ann.mix(VCF_ANNOTATE_SNPEFF.out.vcf_tbi)
        versions = versions.mix(VCF_ANNOTATE_SNPEFF.out.versions)
    }

    if (tools.split(',').contains('merge')) {
        vcf_ann_for_merge = VCF_ANNOTATE_SNPEFF.out.vcf_tbi.map { meta, vcf_, _tbi -> [meta, vcf_, []] }
        VCF_ANNOTATE_MERGE(vcf_ann_for_merge, fasta, vep_genome, vep_species, vep_cache_version, vep_cache, vep_extra_files)

        reports = reports.mix(VCF_ANNOTATE_MERGE.out.reports)
        vcf_ann = vcf_ann.mix(VCF_ANNOTATE_MERGE.out.vcf_tbi)
        versions = versions.mix(VCF_ANNOTATE_MERGE.out.versions)
    }

    if (tools.split(',').contains('vep')) {
        vcf_for_vep = vcf.map { meta, vcf_ -> [meta, vcf_, []] }
        VCF_ANNOTATE_ENSEMBLVEP(vcf_for_vep, fasta, vep_genome, vep_species, vep_cache_version, vep_cache, vep_extra_files)

        reports = reports.mix(VCF_ANNOTATE_ENSEMBLVEP.out.reports)
        vcf_ann = vcf_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.vcf_tbi)
        tab_ann = tab_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.tab)
        json_ann = json_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.json)
        versions = versions.mix(VCF_ANNOTATE_ENSEMBLVEP.out.versions)
    }

    // SnpSift runs LAST on annotated outputs (per best practices)
    // If no other annotators were used, fall back to original vcf
    if (tools.split(',').contains('snpsift')) {
        def has_other_annotators = tools.split(',').any { it in ['bcfann', 'snpeff', 'vep', 'merge'] }

        if (has_other_annotators) {
            vcf_for_snpsift = vcf_ann.map { meta, vcf_, _tbi -> [meta, vcf_] }
        } else {
            vcf_for_snpsift = vcf
        }

        VCF_ANNOTATE_SNPSIFT(
            vcf_for_snpsift,
            snpsift_db_configs
        )

        vcf_ann = vcf_ann.mix(VCF_ANNOTATE_SNPSIFT.out.vcf_tbi)
        // versions collected via topic channel
    }

    emit:
    vcf_ann  // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    tab_ann
    json_ann
    reports  //    path: *.html
    versions //    path: versions.yml
}
