#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
kate: syntax groovy; space-indent on; indent-width 2;
================================================================================
=                                 S  A  R  E  K                                =
================================================================================
 New Germline (+ Somatic) Analysis Workflow. Started March 2016.
--------------------------------------------------------------------------------
 @Authors
 Sebastian DiLorenzo <sebastian.dilorenzo@bils.se> [@Sebastian-D]
 Jesper Eisfeldt <jesper.eisfeldt@scilifelab.se> [@J35P312]
 Phil Ewels <phil.ewels@scilifelab.se> [@ewels]
 Maxime Garcia <maxime.garcia@scilifelab.se> [@MaxUlysse]
 Szilveszter Juhos <szilveszter.juhos@scilifelab.se> [@szilvajuhos]
 Max Käller <max.kaller@scilifelab.se> [@gulfshores]
 Malin Larsson <malin.larsson@scilifelab.se> [@malinlarsson]
 Marcel Martin <marcel.martin@scilifelab.se> [@marcelm]
 Björn Nystedt <bjorn.nystedt@scilifelab.se> [@bjornnystedt]
 Pall Olason <pall.olason@scilifelab.se> [@pallolason]
 Pelin Sahlén <pelin.akan@scilifelab.se> [@pelinakan]
--------------------------------------------------------------------------------
 @Homepage
 http://opensource.scilifelab.se/projects/sarek/
--------------------------------------------------------------------------------
 @Documentation
 https://github.com/SciLifeLab/Sarek/README.md
--------------------------------------------------------------------------------
 Processes overview
 - GenerateMultiQCconfig - Generate MultiQC configuration file
 - RunMultiQC - Run MultiQC on reports
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
try {
    if( ! nextflow.version.matches(">= ${params.nfRequiredVersion}") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version ${params.nfRequiredVersion} required! You are running v${workflow.nextflow.version}.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please update Nextflow.\n" +
              "============================================================"
}

if (params.help) exit 0, helpMessage()
if (!SarekUtils.isAllowedParams(params)) exit 1, "params unknown, see --help for more information"
if (!checkUppmaxProject()) exit 1, "No UPPMAX project ID found! Use --project <UPPMAX Project ID>"

directoryMap = defineDirectoryMap()
/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

startMessage()

process GenerateMultiQCconfig {
  publishDir directoryMap.multiQC, mode: 'link'

  input:

  output:
  file("multiqc_config.yaml") into multiQCconfig

  when: !params.noReports

  script:
  """
  touch multiqc_config.yaml
  echo "custom_logo: ${baseDir}/doc/images/Sarek_no_Border.png" >> multiqc_config.yaml
  echo "custom_logo_url: http://opensource.scilifelab.se/projects/sarek" >> multiqc_config.yaml
  echo "custom_logo_title: 'Sarek'" >> multiqc_config.yaml
  echo "report_header_info:" >> multiqc_config.yaml
  echo "- Sarek version: ${params.version}" >> multiqc_config.yaml
  echo "- Contact Name: ${params.callName}" >> multiqc_config.yaml
  echo "- Contact E-mail: ${params.contactMail}" >> multiqc_config.yaml
  echo "- Directory: ${workflow.launchDir}" >> multiqc_config.yaml
  echo "- Genome: "${params.genome} >> multiqc_config.yaml
  echo "top_modules:" >> multiqc_config.yaml
  echo "- 'fastqc'" >> multiqc_config.yaml
  echo "- 'picard'" >> multiqc_config.yaml
  echo "- 'samtools'" >> multiqc_config.yaml
  echo "- 'qualimap'" >> multiqc_config.yaml
  echo "- 'snpeff'" >> multiqc_config.yaml
  echo "- 'vep'" >> multiqc_config.yaml
  """
}

if (params.verbose && !params.noReports) multiQCconfig = multiQCconfig.view {
  "MultiQC config:\n\
  File  : [${it.fileName}]"
}

reportsForMultiQC = Channel.empty()
  .mix(
    Channel.fromPath("${directoryMap.bamQC}/*", type: 'dir'),
    Channel.fromPath("${directoryMap.bcftoolsStats}/*"),
    Channel.fromPath("${directoryMap.fastQC}/*/*"),
    Channel.fromPath("${directoryMap.markDuplicatesQC}/*"),
    Channel.fromPath("${directoryMap.samtoolsStats}/*"),
    multiQCconfig
  ).collect()

process RunMultiQC {
  publishDir directoryMap.multiQC, mode: 'link'

  input:
    file ('*') from reportsForMultiQC

  output:
    set file("*multiqc_report.html"), file("*multiqc_data") into multiQCReport

  when: !params.noReports

  script:
  """
  multiqc -f -v .
  """
}

if (params.verbose) multiQCReport = multiQCReport.view {
  "MultiQC report:\n\
  File  : [${it[0].fileName}]\n\
  Dir   : [${it[1].fileName}]"
}

/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def checkUppmaxProject() {
  // check if UPPMAX project number is specified
  return !(workflow.profile == 'slurm' && !params.project)
}

def defineDirectoryMap() {
  return [
    'bamQC'            : "${params.outDir}/Reports/bamQC",
    'bcftoolsStats'    : "${params.outDir}/Reports/BCFToolsStats",
    'fastQC'           : "${params.outDir}/Reports/FastQC",
    'markDuplicatesQC' : "${params.outDir}/Reports/MarkDuplicates",
    'multiQC'          : "${params.outDir}/Reports/MultiQC",
    'samtoolsStats'    : "${params.outDir}/Reports/SamToolsStats"
  ]
}

def grabRevision() {
  // Return the same string executed from github or not
  return workflow.revision ?: workflow.commitId ?: workflow.scriptId.substring(0,10)
}

def helpMessage() {
  // Display help message
  this.sarekMessage()
  log.info "    Usage:"
  log.info "       nextflow run runMultiQC.nf"
  log.info "    --help"
  log.info "       you're reading it"
  log.info "    --verbose"
  log.info "       Adds more verbosity to workflow"
}

def minimalInformationMessage() {
  // Minimal information message
  log.info "Command Line: " + workflow.commandLine
  log.info "Profile     : " + workflow.profile
  log.info "Project Dir : " + workflow.projectDir
  log.info "Launch Dir  : " + workflow.launchDir
  log.info "Work Dir    : " + workflow.workDir
  log.info "Out Dir     : " + params.outDir
  log.info "Genome      : " + params.genome
  log.info "Genome_base : " + params.genome_base
  log.info "Containers  :"
  if (params.repository) log.info "  Repository   : ${params.repository}"
  else log.info "  ContainerPath: " + params.containerPath
  log.info "  Tag          : " + params.tag
}

def nextflowMessage() {
  // Nextflow message (version + build)
  log.info "N E X T F L O W  ~  version ${workflow.nextflow.version} ${workflow.nextflow.build}"
}

def sarekMessage() {
  // Display Sarek message
  log.info "Sarek - Workflow For Somatic And Germline Variations ~ ${params.version} - " + this.grabRevision() + (workflow.commitId ? " [${workflow.commitId}]" : "")
}

def startMessage() {
  // Display start message
  this.sarekMessage()
  this.minimalInformationMessage()
}

workflow.onComplete {
  // Display complete message
  this.nextflowMessage()
  this.sarekMessage()
  this.minimalInformationMessage()
  log.info "Completed at: " + workflow.complete
  log.info "Duration    : " + workflow.duration
  log.info "Success     : " + workflow.success
  log.info "Exit status : " + workflow.exitStatus
  log.info "Error report: " + (workflow.errorReport ?: '-')
}

workflow.onError {
  // Display error message
  this.nextflowMessage()
  this.sarekMessage()
  log.info "Workflow execution stopped with the following message:"
  log.info "  " + workflow.errorMessage
}
