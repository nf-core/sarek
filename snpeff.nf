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
<<<<<<< HEAD
  -jar ${params.snpeffHome}/snpEff.jar \
  ${params.snpeffDb} \
=======
  -jar snpEff.jar \
  ${params.snpeffHome} \
>>>>>>> c2ad2ee852cae737bd34eaa11102415fb57d7cdc
  ${params.sampleVcf} > snpeffout.ann.vcf
  """
}
