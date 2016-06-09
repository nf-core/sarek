#!/usr/bin/env nextflow

/*
 * Sample run data for Strelka
 * use like (on milou):
 * ./nextflow run Strelka.nf --tumor_bam ~/dev/chr17_testdata/HCC1143.tumor.bam --normal_bam ~/dev/chr17_testdata/HCC1143.normal.bam
 */

tumor_bam = file(params.tumor_bam)
tumor_bai = file(params.tumor_bai)
normal_bam = file(params.normal_bam)
normal_bai = file(params.normal_bai)
genome_file = file(params.genome)
genome_index = file(params.genomeidx)

process Strelka {
    module 'bioinfo-tools'

    cpus 1

    input:
    file genome_file
    file genome_index
    file tumor_bam
    file tumor_bai
    file normal_bam
    file normal_bai


    """
    """
}
