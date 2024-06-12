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
    bcfann_enabled           // boolean
    bcftools_annotations     // channel (value)
    bcftools_annotations_tbi // channel (value)
    bcftools_header_lines    // channel (value)
    snpeff_enabled           // boolean
    snpeff_db                // resolved value
    snpeff_cache             // channel (value):
    merge_enabled            // boolean
    vep_enabled              // boolean
    vep_cache                // channel (value):
    vep_cache_version        // params
    vep_extra_files          // array?
    vep_fasta                // channel (value)
    vep_genome               // params
    vep_species              // params

    main:
    reports  = Channel.empty()
    vcf_ann  = Channel.empty()
    tab_ann  = Channel.empty()
    json_ann = Channel.empty()
    versions = Channel.empty()

    if (bcfann_enabled) {
        VCF_ANNOTATE_BCFTOOLS(
            input,
            bcftools_annotations,
            bcftools_annotations_tbi,
            bcftools_header_lines)

        vcf_ann  = vcf_ann.mix(VCF_ANNOTATE_BCFTOOLS.out.vcf_tbi)
        versions = versions.mix(VCF_ANNOTATE_BCFTOOLS.out.versions)
    }

    if (snpeff_enabled) {
        VCF_ANNOTATE_SNPEFF(
            input,
            snpeff_db,
            snpeff_cache)

        reports  = reports.mix(VCF_ANNOTATE_SNPEFF.out.reports.map{ meta, reports -> [ reports ] })
        vcf_ann  = vcf_ann.mix(VCF_ANNOTATE_SNPEFF.out.vcf_tbi)
        versions = versions.mix(VCF_ANNOTATE_SNPEFF.out.versions)
    }

    if (merge_enabled) {
        VCF_ANNOTATE_MERGE(
            VCF_ANNOTATE_SNPEFF.out.vcf_tbi.map{ meta, vcf, tbi -> [ meta, vcf, [] ] },
            vep_fasta,
            vep_genome,
            vep_species,
            vep_cache_version,
            vep_cache,
            vep_extra_files)

        reports  = reports.mix(VCF_ANNOTATE_MERGE.out.reports)
        vcf_ann  = vcf_ann.mix(VCF_ANNOTATE_MERGE.out.vcf_tbi)
        versions = versions.mix(VCF_ANNOTATE_MERGE.out.versions)
    }

    if (vep_enabled) {
        VCF_ANNOTATE_ENSEMBLVEP(
            input.map{ meta, vcf -> [ meta, vcf, [] ] },
            vep_fasta,
            vep_genome,
            vep_species,
            vep_cache_version,
            vep_cache,
            vep_extra_files)

        reports  = reports.mix(VCF_ANNOTATE_ENSEMBLVEP.out.reports)
        vcf_ann  = vcf_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.vcf_tbi)
        tab_ann  = tab_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.tab)
        json_ann = json_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.json)
        versions = versions.mix(VCF_ANNOTATE_ENSEMBLVEP.out.versions)
    }

    emit:
    vcf      = vcf_ann  // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    tab      = tab_ann  // only from VEP
    json     = json_ann // only from VEP
    reports             //    path: *.html
    versions            //    path: versions.yml
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
