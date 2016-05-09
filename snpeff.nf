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
  -jar snpEff.jar \
  ${params.snpeffHome} \
  ${params.sampleVcf} > snpeffout.ann.vcf
  """
}
