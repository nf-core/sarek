#!/usr/bin/env nextflow

/*
 * Sample run data for Variant Calling
 * use like (on milou):
 * interactive -p devel -A b2013064 -t 60 nextflow -c milou.config run variantCalling.nf --tumorBam ~/dev/chr17_testdata/HCC1143.tumor.bam --normalBam ~/dev/chr17_testdata/HCC1143.normal.bam
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
  input:
  file tumorBam
  file normalBam

  output:
  file mutect1Vcf

  """
  nextflow run $baseDir/Mutect1.nf \
  -w ${params.cawDir} \
  --tumorBam ${params.tumorBam} \
  --normalBam ${params.normalBam} \
  --genome ${params.genome} \
  --cosmic ${params.cosmic} \
  --dbsnp ${params.dbsnp} \
  --mutect1Home ${params.mutect1Home}
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
  input:
  file tumorBam
  file normalBam

  output:
  file vardictVcf

  """
  nextflow run $baseDir/Vardict.nf \
  -w ${params.cawDir} \
  --tumorBam ${params.tumorBam} \
  --normalBam ${params.normalBam} \
  --genome ${params.genome} \
  --cosmic ${params.cosmic} \
  --dbsnp ${params.dbsnp}
  """
}