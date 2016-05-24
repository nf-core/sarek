#!/usr/bin/env nextflow

/*
 * Sample run data for VarDictJava
 * use like (on milou):
 * ./nextflow run Vardict.nf --tumorBam ~/dev/chr17_testdata/HCC1143.tumor.bam --normalBam ~/dev/chr17_testdata/HCC1143.normal.bam
 */

tumorBam = file(params.tumorBam)
tumorBai = file(params.tumorBai)
normalBam = file(params.normalBam)
normalBai = file(params.normalBai)
genomeFile = file(params.genome)
genomeIndex = file(params.genomeidx)

process Vardict {
  module 'bioinfo-tools'
  module 'VarDictJava/1.4.5'

  cpus 1

  input:
  file genomeFile
  file genomeIndex

  file tumorBam
  file tumorBai
  file normalBam
  file normalBai

  output:
  file '*.VarDict.vcf' into vardict_vc_call

  """
  VarDict -G ${params.genome} \
  -f 0.01 -N TSN \
  -b "${normalBam}|${tumorBam}" \
  -z 1 -F 0x500 \
  -c 1 -S 2 -E 3 -g 4 \
  -R chr17:1000000-1100000 | testsomatic.R | var2vcf_somatic.pl > test.VarDict.vcf
  """
}