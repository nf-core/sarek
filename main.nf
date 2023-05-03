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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.ascat_alleles         = WorkflowMain.getGenomeAttribute(params, 'ascat_alleles')
params.ascat_genome          = WorkflowMain.getGenomeAttribute(params, 'ascat_genome')
params.ascat_loci            = WorkflowMain.getGenomeAttribute(params, 'ascat_loci')
params.ascat_loci_gc         = WorkflowMain.getGenomeAttribute(params, 'ascat_loci_gc')
params.ascat_loci_rt         = WorkflowMain.getGenomeAttribute(params, 'ascat_loci_rt')
params.bwa                   = WorkflowMain.getGenomeAttribute(params, 'bwa')
params.bwamem2               = WorkflowMain.getGenomeAttribute(params, 'bwamem2')
params.cf_chrom_len          = WorkflowMain.getGenomeAttribute(params, 'cf_chrom_len')
params.chr_dir               = WorkflowMain.getGenomeAttribute(params, 'chr_dir')
params.dbsnp                 = WorkflowMain.getGenomeAttribute(params, 'dbsnp')
params.dbsnp_tbi             = WorkflowMain.getGenomeAttribute(params, 'dbsnp_tbi')
params.dbsnp_vqsr            = WorkflowMain.getGenomeAttribute(params, 'dbsnp_vqsr')
params.dict                  = WorkflowMain.getGenomeAttribute(params, 'dict')
params.dragmap               = WorkflowMain.getGenomeAttribute(params, 'dragmap')
params.fasta                 = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.fasta_fai             = WorkflowMain.getGenomeAttribute(params, 'fasta_fai')
params.germline_resource     = WorkflowMain.getGenomeAttribute(params, 'germline_resource')
params.germline_resource_tbi = WorkflowMain.getGenomeAttribute(params, 'germline_resource_tbi')
params.intervals             = WorkflowMain.getGenomeAttribute(params, 'intervals')
params.known_snps            = WorkflowMain.getGenomeAttribute(params, 'known_snps')
params.known_snps_tbi        = WorkflowMain.getGenomeAttribute(params, 'known_snps_tbi')
params.known_snps_vqsr       = WorkflowMain.getGenomeAttribute(params, 'known_snps_vqsr')
params.known_indels          = WorkflowMain.getGenomeAttribute(params, 'known_indels')
params.known_indels_tbi      = WorkflowMain.getGenomeAttribute(params, 'known_indels_tbi')
params.known_indels_vqsr     = WorkflowMain.getGenomeAttribute(params, 'known_indels_vqsr')
params.mappability           = WorkflowMain.getGenomeAttribute(params, 'mappability')
params.pon                   = WorkflowMain.getGenomeAttribute(params, 'pon')
params.pon_tbi               = WorkflowMain.getGenomeAttribute(params, 'pon_tbi')
params.snpeff_db             = WorkflowMain.getGenomeAttribute(params, 'snpeff_db')
params.snpeff_genome         = WorkflowMain.getGenomeAttribute(params, 'snpeff_genome')
params.snpeff_version        = WorkflowMain.getGenomeAttribute(params, 'snpeff_version')
params.vep_cache_version     = WorkflowMain.getGenomeAttribute(params, 'vep_cache_version')
params.vep_genome            = WorkflowMain.getGenomeAttribute(params, 'vep_genome')
params.vep_species           = WorkflowMain.getGenomeAttribute(params, 'vep_species')
params.vep_version           = WorkflowMain.getGenomeAttribute(params, 'vep_version')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAREK } from './workflows/sarek'

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
    NFCORE_SAREK ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
