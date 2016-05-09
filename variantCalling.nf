#!/usr/bin/env nextflow

/*
 * Sample run data for Variant Calling
 * use like (on milou):
 * ./nextflow run variantCalling.nf --tumor_bam ~/dev/chr17_testdata/HCC1143.tumor.bam --normal_bam ~/dev/chr17_testdata/HCC1143.normal.bam
 */

version="0.0.1"

// First checking the existence of the BAM file and its index for the tumor
params.tumor_bam = "tumor.bam"
tumor_bam = file(params.tumor_bam)
if(!tumor_bam.exists()) exit 1, "Missing tumor file ${tumor_bam}; please specify --tumor_bam <SAMPLE> --normal_bam <SAMPLE> " 
// the index
params.tumor_bai = params.tumor_bam.replaceFirst(/.bam/,".bam.bai")
tumor_bai = file(params.tumor_bai)
if(!tumor_bai.exists()) exit 1, "Missing tumor file index ${tumor_bai}; please run samtools index ${normal_bam}" 

// Ditto for the normal
params.normal_bam = "normal.bam"
normal_bam = file(params.normal_bam)
if(!normal_bam.exists()) exit 1, "Missing normal file ${normal_bam}; please specify --tumor_bam <SAMPLE> --normal_bam <SAMPLE>"  
// the normal index
params.normal_bai = params.normal_bam.replaceFirst(/.bam/,".bam.bai")
normal_bai = file(params.normal_bai)
if(!normal_bai.exists()) exit 1, "Missing normal file index ${normal_bai}; please run samtools index ${normal_bam}"

genomeFile  = file(params.genome)
genomeIndex = file(params.genomeIndex)
cosmic      = file(params.cosmic)
cosmicIndex = file(params.cosmicIndex)
dbsnp       = file(params.dbsnp)
dbsnpIndex  = file(params.dbsnpIndex)

// Basic argument handling

switch (params){
  case help :
    text = Channel.from(
      "CANCER ANALYSIS WORKFLOW ~ version $version",
      "    --help",
      "       you're reading it",
      "    --version",
      "       displays version number")
    text.subscribe { println "$it" }
    exit 1

  case version :
    text = Channel.from(
      "CANCER ANALYSIS WORKFLOW",
      "  Version $version",
      "  Last update on $dateUpdate",
      "Project : $workflow.projectDir",
      "Git info: $workflow.repository - $workflow.revision [$workflow.commitId]",
      "Cmd line: $workflow.commandLine")
    text.subscribe { println "$it" }
    exit 1
}

process mutect1 {
  """
  nextflow run $baseDir/MuTect1.nf \
  -w ${params.cawDir} \
  --tumor_bam ${params.tumor_bam} \
  --normal_bam ${params.normal_bam} \
  --genome ${params.genome} \
  --genomeidx ${params.genomeIndex} \
  --cosmic ${params.cosmic} \
  --cosmicidx ${params.cosmicIndex} \
  --dbsnp ${params.dbsnp} \
  --dbsnpidx ${params.dbsnpIndex} \
  --mutect1_home ${params.mutect1Home}
  """
}

// process FreeBayes {
//  """
//  nextflow run $baseDir/FreeBayes.nf \
//  -w ~/workspace/nextflow/work \
//  --tumor_bam params.tumor_bam
//  --normal_bam params.normal_bam
//  """
// }

// process VarDictJava {
//  """
//  nextflow run $baseDir/VarDictJava.nf \
//  -w ~/workspace/nextflow/work \
//  --tumor_bam params.tumor_bam
//  --normal_bam params.normal_bam
//  """
// }

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