#!/usr/bin/env nextflow

/*
--------------------------------------------------------------------------------
                                  nf-core/sarek
--------------------------------------------------------------------------------
Started March 2016.
Ported to nf-core May 2019.
Ported to DSL 2 July 2020.
--------------------------------------------------------------------------------
nf-core/sarek:
  An open-source analysis pipeline to detect germline or somatic variants
  from whole genome or targeted sequencing
--------------------------------------------------------------------------------
 @Homepage
 https://nf-co.re/sarek
--------------------------------------------------------------------------------
 @Documentation
 https://nf-co.re/sarek/latest/usage
--------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////

def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/sarek -profile docker --input sample.tsv --genome GRCh38"
    log.info Schema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --          PARAMETER CHECKS                -- */
////////////////////////////////////////////////////

// Check that conda channels are set-up correctly
if (params.enable_conda) {
    Checks.check_conda_channels(log)
}

// Check AWS batch settings
Checks.aws_batch(workflow, params)

// Check the hostnames against configured profiles
Checks.hostname(workflow, params, log)

// Check genome key exists if provided
Checks.genome_exists(params, log)

////////////////////////////////////////////////////
/* --        REFERENCES PARAMETER VALUES       -- */
////////////////////////////////////////////////////
/* -- Initialize each params in params.genomes -- */
/* --  catch the command line first if defined -- */
////////////////////////////////////////////////////

params.ac_loci                 = Checks.get_genome_attribute(params, 'ac_loci')
params.ac_loci_gc              = Checks.get_genome_attribute(params, 'ac_loci_gc')
params.bwa                     = Checks.get_genome_attribute(params, 'bwa')
params.chr_dir                 = Checks.get_genome_attribute(params, 'chr_dir')
params.chr_length              = Checks.get_genome_attribute(params, 'chr_length')
params.dbsnp                   = Checks.get_genome_attribute(params, 'dbsnp')
params.dbsnp_index             = Checks.get_genome_attribute(params, 'dbsnp_index')
params.dict                    = Checks.get_genome_attribute(params, 'dict')
params.fasta                   = Checks.get_genome_attribute(params, 'fasta')
params.fasta_fai               = Checks.get_genome_attribute(params, 'fasta_fai')
params.germline_resource       = Checks.get_genome_attribute(params, 'germline_resource')
params.germline_resource_index = Checks.get_genome_attribute(params, 'germline_resource_index')
params.intervals               = Checks.get_genome_attribute(params, 'intervals')
params.known_indels            = Checks.get_genome_attribute(params, 'known_indels')
params.known_indels_index      = Checks.get_genome_attribute(params, 'known_indels_index')
params.mappability             = Checks.get_genome_attribute(params, 'mappability')
params.snpeff_db               = Checks.get_genome_attribute(params, 'snpeff_db')
params.species                 = Checks.get_genome_attribute(params, 'species')
params.vep_cache_version       = Checks.get_genome_attribute(params, 'vep_cache_version')

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////

def summary_params = Schema.params_summary_map(workflow, params, json_schema)
log.info Schema.params_summary_log(workflow, params, json_schema)

////////////////////////////////////////////////////
/* --             RUN THE WORKFLOW             -- */
////////////////////////////////////////////////////

workflow {

    include { SAREK } from './workflows/sarek' addParams( summary_params: summary_params )
    SAREK ()

}
