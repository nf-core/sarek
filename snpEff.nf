#!/usr/bin/env nextflow
/*
 * Sample run data for snpEff
 * use like (on milou):
 * nextflow run snpEff.nf --vcf my.vcf
 */

// checking the existence of the input VCF 
params.vcf= "test.vcf" // override with --vcf <infile>
vcf_call = file(params.vcf)
if(!vcf_call.exists()) exit 1, "Missing variant call results file ${vcf_call}" 

params.genome_version = "GRCh37.75"     // in theory you can change it (using --genome_version <VERSION>) if you want to  

process snpEff_annotation {
    module 'bioinfo-tools'
    module 'snpEff/4.2'

    cpus 1

    input:
    file vcf_call

    output:
    file '*.snpEff.vcf' into snpEff_annotated

    """
    snpEff ${params.genome_version} ${vcf_call} > annotated.snpEff.vcf
    """
}
