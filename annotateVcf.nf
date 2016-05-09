#!/usr/bin/env nextflow

/*
 * Sample run data for Annotation of results
 * use like (on milou):
 * ./nextflow run annotateVcf.nf
 */

version="0.0.1"

snpeffHome = "${params.snpeffHome}"
snpeffDb   = file(params.snpeffDb)

// Basic argument handling

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
      "Project : $workflow.projectDir",
      "Git info: $workflow.repository - $workflow.revision [$workflow.commitId]",
      "Cmd line: $workflow.commandLine")
    text.subscribe { println "$it" }
    exit 1
}

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

process snpeff {
  """
  nextflow run $baseDir/snpeff.nf \
  -w ${params.cawDir} \
  --snpeffHome ${params.snpeffHome}
  --snpeffDb ${params.snpeffDb}
  --sampleVcf ${params.sampleVcf}
  """
}

// process vcfanno {
//   """
//   nextflow run $baseDir/vcfanno.nf \
//   -w ${params.cawDir} \

//   """
// }