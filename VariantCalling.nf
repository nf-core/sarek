#!/usr/bin/env nextflow

/*
 * Sample run data for Variant Calling
 * use like (on milou):
 * interactive -p devel -A b2013064 -t 60 nextflow -c milou.config run VariantCalling.nf --tumorBam ~/dev/data/HCC1143.tumor.bam --normalBam ~/dev/data/HCC1143.normal.bam
 */

String version="0.0.1"
String dateUpdate = "2016-05-24"

/*
 * Get some basic informations about the workflow
 * to get more informations use -with-trace or -with-timeline
 */

workflow.onComplete {
  text = Channel.from(
    "CANCER ANALYSIS WORKFLOW",
    "Version     : $version",
    "Command line: ${workflow.commandLine}",
    "Completed at: ${workflow.complete}",
    "Duration    : ${workflow.duration}",
    "Success     : ${workflow.success}",
    "workDir     : ${workflow.workDir}",
    "Exit status : ${workflow.exitStatus}",
    "Error report: ${workflow.errorReport ?: '-'}")
  text.subscribe { log.info "$it" }
}

/*
 * Basic argument handling
 * and verification
 */

switch (params) {
  case {params.help} :
    text = Channel.from(
      "CANCER ANALYSIS WORKFLOW ~ version $version",
      "    --help",
      "       you're reading it",
      "    --version",
      "       displays version number")
    text.subscribe { println "$it" }
    exit 1

  case {params.version} :
    text = Channel.from(
      "CANCER ANALYSIS WORKFLOW",
      "  Version $version",
      "  Last update on $dateUpdate",
      "Project : $workflow.projectDir",
      "Cmd line: $workflow.commandLine")
    text.subscribe { println "$it" }
    exit 1

  case {params.tumorBam} :
    tumorBam = file(params.tumorBam)
    params.tumorBai = params.tumorBam.replaceFirst(/.bam/,".bam.bai")
    tumorBai = file(params.tumorBai)

  case {params.normalBam} :
    normalBam = file(params.normalBam)
    params.normalBai = params.normalBam.replaceFirst(/.bam/,".bam.bai")
    normalBai = file(params.normalBai)

  case {params.genome} :
   genomeFile = file(params.genome)

  case {params.genomeIndex} :
    genomeIndex = file(params.genomeIndex)

  case {params.cosmic} :
   cosmic = file(params.cosmic)

  case {params.cosmicIndex} :
    cosmicIndex = file(params.cosmicIndex)

  case {params.dbsnp} :
   dbsnp = file(params.dbsnp)

  case {params.dbsnpIndex} :
    dbsnpIndex = file(params.dbsnpIndex)  
}

// First checking the existence of the BAM file and its index for the tumor

if(!tumorBam.exists())  exit 1, "Missing tumor file ${tumorBam}; please specify --tumorBam <SAMPLE> --normalBam <SAMPLE> " 
if(!tumorBai.exists())  exit 1, "Missing tumor file index ${tumorBai}; please run samtools index ${normalBam}" 
if(!normalBam.exists()) exit 1, "Missing normal file ${normalBam}; please specify --tumorBam <SAMPLE> --normalBam <SAMPLE>"  
if(!normalBai.exists()) exit 1, "Missing normal file index ${normalBai}; please run samtools index ${normalBam}"

process Mutect1 {

  module 'bioinfo-tools'
  moduule 'mutect/1.1.5'

  input:
  file tumorBam
  file normalBam

  output:
  file '*.mutect1.vcf' into mutect1Vcf
  file '*.mutect1.out' into mutect1Out

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

// process Strelka {
//   input:
//   file tumorBam
//   file normalBam

//   output:
//   file strelkaVcf

//   """
//   nextflow run $baseDir/Strelka.nf \
//   -w ${params.cawDir}
//   """
// }

process Vardict {

  module 'bioinfo-tools'
  module 'VarDictJava/1.4.5'
  module 'samtools/1.3'

  cpus 1

  input:
  file tumorBam
  file normalBam

  output:
  file '*.VarDict.vcf' into vardictVcf


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
  ${params.vardictHome}/testsomatic.R | \
  ${params.vardictHome}/var2vcf_somatic.pl -f 0.01 -N "${tumorBam}|${normalBam}" > test.VarDict.vcf
  """
}

process MergeVcf {

  input:
  file mutect1Vcf
  file vardictVcf

  """
  java -Xmx7g
  -jar ${picardHome}/MergeVcfs.jar \
  ${mutect1Vcf} \
  ${vardictVcf}
  """
}