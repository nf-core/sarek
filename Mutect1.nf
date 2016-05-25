#!/usr/bin/env nextflow

/*
 * Sample run data for MuTect1
 * use like (on milou):
 * ./nextflow run Mutect1.nf --tumorBam ~/dev/chr17_testdata/HCC1143.tumor.bam --normalBam ~/dev/chr17_testdata/HCC1143.normal.bam
 */

tumorBam    = file(params.tumorBam)
normalBam   = file(params.normalBam)
genomeFile  = file(params.genome)
genomeIndex = file(params.genomeIndex)
cosmic      = file(params.cosmic)
cosmicIndex = file(params.cosmicIndex)
dbsnp       = file(params.dbsnp)
dbsnpIndex  = file(params.dbsnpIndex)
mutect1Home = ${params.mutect1Home}

process Mutect1 {

  cpus 2

  input:
  file genomeFile
  file cosmic
  file dbsnp
  file genomeIndex
  file tumorBam
  file normalBam

  output:
  file '*.mutect1.vcf' into mutect1Vcf
  file '*.mutect1.out' into mutect1Out

  """
  java -jar ${params.mutect1Home}/muTect-1.1.5.jar \
  --analysis_type MuTect \
  --reference_sequence ${genomeFile} \
  --cosmic ${cosmic} \
  --dbsnp ${dbsnp} \
  --input_file:normal ${normalBam} \
  --input_file:tumor ${tumorBam} \
  --out test.mutect1.out \
  --vcf test.mutect1.vcf \
  -L 17:1000000-2000000
  """
}
