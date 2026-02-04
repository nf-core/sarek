//
// ANNOTATION
//

include { BCFTOOLS_ANNOTATE                             } from '../../../modules/nf-core/bcftools/annotate'
include { ENSEMBLVEP_VEP                                } from '../../../modules/nf-core/ensemblvep/vep'
include { ENSEMBLVEP_VEP as VCF_ANNOTATE_MERGE          } from '../../../modules/nf-core/ensemblvep/vep'
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
    snpsift_db                  // channel: [[databases], [tbis], [vardbs], [fields], [prefixes]]

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
        VCF_ANNOTATE_MERGE(vcf_ann_for_merge, vep_genome,vep_species,vep_cache_version, vep_cache, fasta, vep_extra_files)

        vcf_ann = vcf_ann.mix(VCF_ANNOTATE_MERGE.out.vcf.join(VCF_ANNOTATE_MERGE.out.tbi, failOnDuplicate: true, failOnMismatch: true))
    }

    if (tools.split(',').contains('vep')) {
        vcf_for_vep = vcf.map { meta, vcf_ -> [meta, vcf_, []] }
        ENSEMBLVEP_VEP(vcf_for_vep, vep_genome, vep_species, vep_cache_version, vep_cache, fasta, vep_extra_files)

        vcf_ann = vcf_ann.mix(ENSEMBLVEP_VEP.out.vcf.join(ENSEMBLVEP_VEP.out.tbi, failOnDuplicate: true, failOnMismatch: true))
        tab_ann = tab_ann.mix(ENSEMBLVEP_VEP.out.tab)
        json_ann = json_ann.mix(ENSEMBLVEP_VEP.out.json)
    }

    // SnpSift runs on all final annotated outputs
    // If no other annotators were used, fall back to original vcf
    if (tools.split(',').contains('snpsift')) {
        def has_other_annotators = ['merge', 'snpeff', 'vep', 'bcfann'].any { tools.split(',').contains(it) }
        def snpsift_input = tools.split(',').contains('merge')
            ? VCF_ANNOTATE_MERGE.out.vcf
            : (has_other_annotators ? vcf_ann : vcf)

        VCF_ANNOTATE_SNPSIFT(snpsift_input, snpsift_db)
        vcf_ann = vcf_ann.mix(VCF_ANNOTATE_SNPSIFT.out.vcf_tbi)
    }

    emit:
    vcf_ann  // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    tab_ann
    json_ann
    reports  //    path: *.html
    versions //    path: versions.yml
}
