#!/usr/bin/env nextflow

/*
 * Sample run data for MuTect1
 * use like (on milou):
 * ./nextflow run Mutect1.nf --tumorBam ~/dev/chr17_testdata/HCC1143.tumor.bam --normalBam ~/dev/chr17_testdata/HCC1143.normal.bam
 */

process Mutect1 {

  cpus 2

  """
  java -jar ${params.mutect1Home}/muTect-1.1.5.jar \
  --analysis_type MuTect \
  --reference_sequence ${params.genome} \
  --cosmic ${params.cosmic} \
  --dbsnp ${params.dbsnp} \
  --input_file:normal ${params.normalBam} \
  --input_file:tumor ${params.tumorBam} \
  --out test.mutect1.out \
  --vcf test.mutect1.vcf \
  -L 17:1000000-2000000
  """
}
