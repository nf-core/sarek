#!/usr/bin/env nextflow

/*
 * Sample run data for VarDictJava
 * use like (on milou):
 * ./nextflow run Vardict.nf --tumorBam ~/dev/chr17_testdata/HCC1143.tumor.bam --normalBam ~/dev/chr17_testdata/HCC1143.normal.bam
 * The original command line to be used should be something like
    AF_THR="0.01" # minimum allele frequency
    <path_to_vardict_folder>/build/install/VarDict/bin/VarDict -G /path/to/hg19.fa \
        -f $AF_THR -N tumor_sample_name -b "/path/to/tumor.bam|/path/to/normal.bam" \
        -z -F -c 1 -S 2 -E 3 -g 4 /path/to/my.bed | VarDict/testsomatic.R | \
        VarDict/var2vcf_somatic.pl -N "tumor_sample_name|normal_sample_name" -f $AF_THR

    In practice:

     /sw/apps/bioinfo/VarDictJava/1.4.5/milou/VarDictJava/build/install/VarDict/bin/VarDict -G /sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/human_g1k_v37_decoy.fasta -f 0.01 -N TSN -b "/home/szilva/glob/chr17_testdata/HCC1143.normal.bam|/home/szilva/glob/chr17_testdata/HCC1143.tumor.bam" -z 1 -F 0x500 -c 1 -S 2 -E 3 -g 4 -R chr17:1000000-1100000 | /sw/apps/bioinfo/VarDictJava/1.4.5/milou/VarDictJava/VarDict/testsomatic.R | /sw/apps/bioinfo/VarDictJava/1.4.5/milou/VarDictJava/VarDict/var2vcf_somatic.pl| tee vdj.vcf

 */

tumorBam  = file (params.tumorBam)
normalBam = file (params.normalBam)

process Vardict {

  module 'bioinfo-tools'
  module 'samtools/1.3'
  module 'VarDictJava/1.4.5'

  cpus 1

  input:
  file tumorBam
  file normalBam

  output:
  file '*.VarDict.vcf' into vardict_vc_call


  // perl /sw/apps/bioinfo/VarDictJava/1.4.5/milou/VarDictJava/vardict.pl -G ${params.genome} \
  // -f 0.01 -N TSN \
  // -b "${params.normalBam}|${params.tumorBam}" \
  // -z 1 -F 0x500 \
  // -c 1 -S 2 -E 3 -g 4 \
  // -R chr17:1000000-1100000 | testsomatic.R | var2vcf_somatic.pl > test.VarDict.vcf

  """
  vardict -G ${params.genome} \
  -f 0.01 -N ${tumorBam} \
  -b "${params.tumorBam}|${params.normalBam}" \
  -z 1 -F 0x500 \
  -c 1 -S 2 -E 3 -g 4 \
  -R chr17:1000000-1100000 | \
  {params.vardictHome}/testsomatic.R | \
  {params.vardictHome}/var2vcf_somatic.pl -f 0.01 -N "${tumorBam}|${normalBam} > test.VarDict.vcf
  """
}
