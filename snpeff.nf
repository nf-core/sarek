#!/usr/bin/env nextflow

/*
 * Sample run data for snpEff
 * use like (on milou):
 * ./nextflow run snpeff.nf
 */

vcfFile = file(params.sampleVcf)

process Snpeff {

  cpus 2

  input:
  file vcfFile

  output:
  file '*.ann.vcf' into snpeffOut

  """
  java -Xmx4g \
  -jar ${params.snpeffHome}/snpEff.jar \
  ${params.snpeffDb} \
  ${params.sampleVcf} > snpeffout.ann.vcf
  """
}
