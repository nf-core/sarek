#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/sarek
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Started March 2016.
    Ported to nf-core May 2019.
    Ported to DSL 2 July 2020.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/sarek:
        An open-source analysis pipeline to detect germline or somatic variants
        from whole genome or targeted sequencing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/sarek
    Website: https://nf-co.re/sarek
    Docs   : https://nf-co.re/sarek/usage
    Slack  : https://nfcore.slack.com/channels/sarek
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

include { getGenomeAttribute } from './subworkflows/local/utils_sarek'
include { retrieveInput      } from './subworkflows/local/utils_sarek'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.ascat_alleles           = getGenomeAttribute(params, 'ascat_alleles')
params.ascat_genome            = getGenomeAttribute(params, 'ascat_genome')
params.ascat_loci              = getGenomeAttribute(params, 'ascat_loci')
params.ascat_loci_gc           = getGenomeAttribute(params, 'ascat_loci_gc')
params.ascat_loci_rt           = getGenomeAttribute(params, 'ascat_loci_rt')
params.bwa                     = getGenomeAttribute(params, 'bwa')
params.bwamem2                 = getGenomeAttribute(params, 'bwamem2')
params.cf_chrom_len            = getGenomeAttribute(params, 'cf_chrom_len')
params.chr_dir                 = getGenomeAttribute(params, 'chr_dir')
params.dbsnp                   = getGenomeAttribute(params, 'dbsnp')
params.dbsnp_tbi               = getGenomeAttribute(params, 'dbsnp_tbi')
params.dbsnp_vqsr              = getGenomeAttribute(params, 'dbsnp_vqsr')
params.dict                    = getGenomeAttribute(params, 'dict')
params.dragmap                 = getGenomeAttribute(params, 'dragmap')
params.fasta                   = getGenomeAttribute(params, 'fasta')
params.fasta_fai               = getGenomeAttribute(params, 'fasta_fai')
params.germline_resource       = getGenomeAttribute(params, 'germline_resource')
params.germline_resource_tbi   = getGenomeAttribute(params, 'germline_resource_tbi')
params.intervals               = getGenomeAttribute(params, 'intervals')
params.known_indels            = getGenomeAttribute(params, 'known_indels')
params.known_indels_tbi        = getGenomeAttribute(params, 'known_indels_tbi')
params.known_indels_vqsr       = getGenomeAttribute(params, 'known_indels_vqsr')
params.known_snps              = getGenomeAttribute(params, 'known_snps')
params.known_snps_tbi          = getGenomeAttribute(params, 'known_snps_tbi')
params.known_snps_vqsr         = getGenomeAttribute(params, 'known_snps_vqsr')
params.mappability             = getGenomeAttribute(params, 'mappability')
params.ngscheckmate_bed        = getGenomeAttribute(params, 'ngscheckmate_bed')
params.pon                     = getGenomeAttribute(params, 'pon')
params.pon_tbi                 = getGenomeAttribute(params, 'pon_tbi')
params.sentieon_dnascope_model = getGenomeAttribute(params, 'sentieon_dnascope_model')
params.snpeff_db               = getGenomeAttribute(params, 'snpeff_db')
params.snpeff_genome           = getGenomeAttribute(params, 'snpeff_genome')
params.vep_cache_version       = getGenomeAttribute(params, 'vep_cache_version')
params.vep_genome              = getGenomeAttribute(params, 'vep_genome')
params.vep_species             = getGenomeAttribute(params, 'vep_species')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ALTERNATIVE INPUT FILE ON RESTART
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.input_restart = retrieveInput(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAREK                   } from './workflows/sarek'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_sarek'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_sarek'

// WORKFLOW: Run main nf-core/sarek analysis pipeline
workflow NFCORE_SAREK {

    SAREK ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
workflow {

    PIPELINE_INITIALISATION()

    NFCORE_SAREK ()

    // PIPELINE_COMPLETION()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
