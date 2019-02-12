#!/usr/bin/env nextflow

/*
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
--------------------------------------------------------------------------------
 @Homepage
 http://opensource.scilifelab.se/projects/sarek/
--------------------------------------------------------------------------------
 @Documentation
 https://github.com/SciLifeLab/Sarek/README.md
--------------------------------------------------------------------------------
 Processes overview
 - GetVersionAll - Get version of tools
 - RunMultiQC - Run MultiQC on reports
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

if (params.help) exit 0, helpMessage()
if (!SarekUtils.isAllowedParams(params)) exit 1, "params unknown, see --help for more information"
if (!checkUppmaxProject()) exit 1, "No UPPMAX project ID found! Use --project <UPPMAX Project ID>"

// Check for awsbatch profile configuration
// make sure queue is defined
if (workflow.profile == 'awsbatch') {
    if(!params.awsqueue) exit 1, "Provide the job queue for aws batch!"
}

directoryMap = SarekUtils.defineDirectoryMap(params.outDir)
/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

startMessage()

process GetVersionAll {
  publishDir directoryMap.multiQC, mode: params.publishDirMode

  input:
    file(versions) from Channel.fromPath("${directoryMap.version}/*").collect().ifEmpty(file ("empty"))

  output:
    file ("tool_versions_mqc.yaml") into versionsForMultiQC

  when: !params.noReports

  script:
  """
  bcftools version > v_bcftools.txt 2>&1 || true
  bwa &> v_bwa.txt 2>&1 || true
  configManta.py --version > v_manta.txt 2>&1 || true
  configureStrelkaGermlineWorkflow.py --version > v_strelka.txt 2>&1 || true
  echo "${workflow.manifest.version}" &> v_sarek.txt 2>&1 || true
  echo "${workflow.nextflow.version}" &> v_nextflow.txt 2>&1 || true
  echo "SNPEFF version"\$(snpEff -h 2>&1) > v_snpeff.txt
  fastqc -v > v_fastqc.txt 2>&1 || true
  freebayes --version > v_freebayes.txt 2>&1 || true
  gatk ApplyBQSR --help 2>&1 | grep Version: > v_gatk.txt 2>&1 || true
  multiqc --version &> v_multiqc.txt 2>&1 || true
  qualimap --version &> v_qualimap.txt 2>&1 || true
  samtools --version &> v_samtools.txt 2>&1 || true
  vcftools --version &> v_vcftools.txt 2>&1 || true
  vep --help > v_vep.txt

  scrape_tool_versions.py &> tool_versions_mqc.yaml
  """
}

if (params.verbose && !params.noReports) versionsForMultiQC = versionsForMultiQC.view {
  "MultiQC tools version:\n\
  File  : [${it.fileName}]"
}

reportsForMultiQC = Channel.empty()
  .mix(
    Channel.fromPath("${directoryMap.bamQC}/*", type: 'dir'),
    Channel.fromPath("${directoryMap.bcftoolsStats}/*"),
    Channel.fromPath("${directoryMap.fastQC}/*/*"),
    Channel.fromPath("${directoryMap.markDuplicatesQC}/*"),
    Channel.fromPath("${directoryMap.samtoolsStats}/*"),
    Channel.fromPath("${directoryMap.snpeffReports}/*"),
    Channel.fromPath("${directoryMap.vcftools}/*"),
  ).collect()

process RunMultiQC {
  publishDir directoryMap.multiQC, mode: params.publishDirMode

  input:
    file (multiqcConfig) from createMultiQCconfig()
    file (reports) from reportsForMultiQC
    file (versions) from versionsForMultiQC

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

def createMultiQCconfig() {
  def file = workDir.resolve('multiqc_config.yaml')
  file.text  = """
  custom_logo: ${baseDir}/docs/images/Sarek_no_Border.png
  custom_logo_url: http://opensource.scilifelab.se/projects/sarek
  custom_logo_title: 'Sarek'
  report_header_info:
  - Contact Name: ${params.callName}
  - Contact E-mail: ${params.contactMail}
  - Genome: ${params.genome}
  top_modules:
  - 'fastqc'
  - 'picard'
  - 'samtools'
  - 'qualimap'
  - 'bcftools'
  - 'vcftools'
  - 'snpeff'
  """.stripIndent()

  return file
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
  log.info "Containers"
  if (params.repository != "") log.info "  Repository   : " + params.repository
  if (params.containerPath != "") log.info "  ContainerPath: " + params.containerPath
  log.info "  Tag          : " + params.tag
}

def nextflowMessage() {
  // Nextflow message (version + build)
  log.info "N E X T F L O W  ~  version ${workflow.nextflow.version} ${workflow.nextflow.build}"
}

def sarekMessage() {
  // Display Sarek message
  log.info "Sarek - Workflow For Somatic And Germline Variations ~ ${workflow.manifest.version} - " + this.grabRevision() + (workflow.commitId ? " [${workflow.commitId}]" : "")
}

def startMessage() {
  // Display start message
  SarekUtils.sarek_ascii()
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
