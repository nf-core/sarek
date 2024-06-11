/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VCF_ANNOTATE_BCFTOOLS                         } from '../../subworkflows/local/vcf_annotate_bcftools'
include { VCF_ANNOTATE_ENSEMBLVEP                       } from '../../subworkflows/nf-core/vcf_annotate_ensemblvep'
include { VCF_ANNOTATE_ENSEMBLVEP as VCF_ANNOTATE_MERGE } from '../../subworkflows/nf-core/vcf_annotate_ensemblvep'
include { VCF_ANNOTATE_SNPEFF                           } from '../../subworkflows/nf-core/vcf_annotate_snpeff'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANTANNOTATION {
    take:
    input                    // channel (queue): samplesheet read in from --input
    tools                    // array
    bcftools_annotations     // channel (value):
    bcftools_annotations_tbi // channel (value):
    bcftools_header_lines    // string
    snpeff_db                // string
    snpeff_cache             // channel (value):
    vep_cache                // channel (value):
    vep_cache_version        // string
    vep_extra_files          // array?
    vep_fasta                // channel (value)
    vep_genome               // string
    vep_species              // string

    main:
    reports  = Channel.empty()
    vcf_ann  = Channel.empty()
    tab_ann  = Channel.empty()
    json_ann = Channel.empty()
    versions = Channel.empty()

    if (tools.split(',').contains('bcfann')) {
        VCF_ANNOTATE_BCFTOOLS(input, bcftools_annotations, bcftools_annotations_index, bcftools_header_lines)

        vcf_ann = vcf_ann.mix(VCF_ANNOTATE_BCFTOOLS.out.vcf_tbi)
        versions = versions.mix(VCF_ANNOTATE_BCFTOOLS.out.versions)
    }


    if (tools.split(',').contains('merge') || tools.split(',').contains('snpeff')) {
        VCF_ANNOTATE_SNPEFF(input, snpeff_db, snpeff_cache)

        reports = reports.mix(VCF_ANNOTATE_SNPEFF.out.reports.map{ meta, reports -> [ reports ] })
        vcf_ann = vcf_ann.mix(VCF_ANNOTATE_SNPEFF.out.vcf_tbi)
        versions = versions.mix(VCF_ANNOTATE_SNPEFF.out.versions)
    }

    if (tools.split(',').contains('merge')) {
        VCF_ANNOTATE_MERGE(VCF_ANNOTATE_SNPEFF.out.vcf_tbi.map{ meta, vcf, tbi -> [ meta, vcf, [] ] }, vep_fasta, vep_genome, vep_species, vep_cache_version, vep_cache, vep_extra_files)

        reports = reports.mix(VCF_ANNOTATE_MERGE.out.reports)
        vcf_ann = vcf_ann.mix(VCF_ANNOTATE_MERGE.out.vcf_tbi)
        versions = versions.mix(VCF_ANNOTATE_MERGE.out.versions)
    }

    if (tools.split(',').contains('vep')) {
        VCF_ANNOTATE_ENSEMBLVEP(input.map{ meta, vcf -> [ meta, vcf, [] ] }, vep_fasta, vep_genome, vep_species, vep_cache_version, vep_cache, vep_extra_files)

        reports  = reports.mix(VCF_ANNOTATE_ENSEMBLVEP.out.reports)
        vcf_ann  = vcf_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.vcf_tbi)
        tab_ann  = tab_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.tab)
        json_ann = json_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.json)
        versions = versions.mix(VCF_ANNOTATE_ENSEMBLVEP.out.versions)
    }


    emit:
    vcf      = vcf_ann  // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    tab      = tab_ann
    json     = json_ann
    reports             //    path: *.html
    versions            //    path: versions.yml
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
