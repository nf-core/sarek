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
 - BuildDockerContainers - Build containers using Docker
 - PullSingularityContainers - Pull Singularity containers from Docker Hub
 - PushDockerContainers - Push containers to Docker Hub
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

version = '1.2.5'

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
nf_required_version = '0.25.0'
try {
    if( ! nextflow.version.matches(">= ${nf_required_version}") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version ${nf_required_version} required! You are running v${workflow.nextflow.version}.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please update Nextflow.\n" +
              "============================================================"
}

if (params.help) exit 0, helpMessage()
if (params.version) exit 0, versionMessage()
if (!isAllowedParams(params)) exit 1, "params unknown, see --help for more information"
if (!checkUppmaxProject()) exit 1, "No UPPMAX project ID found! Use --project <UPPMAX Project ID>"

// Default params:
// Such params are overridden by command line or configuration definitions

// containerPath is current Directory
params.containerPath = "${baseDir}"
// all containers to be build
params.containers = 'all'
// Docker will not be used
params.docker = false
// Containers will not be pushed on DockerHub
params.push = false
// DockerHub repository is maxulysse
// TODO Change to a SciLifeLab repository
params.repository = 'maxulysse'
// Singularity will not be used
params.singularity = false

// Define containers to handle (build/push or pull)
containersList = defineContainersList()
containers = params.containers.split(',').collect {it.trim()}
containers = containers == ['all'] ? containersList : containers

// push only to DockerHub, so only when using Docker
push = params.docker && params.push ? true : false

// by default the tag will be the current version
tag = params.tag ? params.tag : version

// to simplify verbose mode
verbose = params.verbose

if (!params.docker && !params.singularity) exit 1, 'No container technology choosed, specify --docker or --singularity, see --help for more information'

if (!checkContainers(containers,containersList)) exit 1, 'Unknown container(s), see --help for more information'

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

startMessage()

dockerContainers = containers
singularityContainers = containers

process BuildDockerContainers {
  tag {"${params.repository}/${container}:${tag}"}

  input:
    val container from dockerContainers

  output:
    val container into containersBuilt

  when: params.docker

  script:
  """
  docker build -t ${params.repository}/${container}:${tag} ${baseDir}/containers/${container}/.
  """
}

if (verbose) containersBuilt = containersBuilt.view {
  "Docker container: ${params.repository}/${it}:${tag} built."
}

process PullSingularityContainers {
  tag {"${params.repository}/${container}:${tag}"}

  publishDir "${params.containerPath}", mode: 'move'

  input:
    val container from singularityContainers

  output:
    file("${container}-${tag}.img") into imagePulled

  when: params.singularity

  script:
  """
  singularity pull --name ${container}-${tag}.img docker://${params.repository}/${container}:${tag}
  """
}

if (verbose) imagePulled = imagePulled.view {
  "Singularity image: ${it.fileName} pulled."
}

process PushDockerContainers {
  tag {params.repository + "/" + container + ":" + tag}

  input:
    val container from containersBuilt

  output:
    val container into containersPushed

  when: params.docker && push

  script:
  """
  docker push ${params.repository}/${container}:${tag}
  """
}

if (verbose) containersPushed = containersPushed.view {
  "Docker container: ${params.repository}/${it}:${tag} pushed."
}

/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def sarekMessage() {
  // Display Sarek message
  log.info "Sarek ~ ${version} - " + this.grabRevision() + (workflow.commitId ? " [${workflow.commitId}]" : "")
}

def checkContainerExistence(container, list) {
  try {assert list.contains(container)}
  catch (AssertionError ae) {
    println("Unknown container: ${container}")
    return false
  }
  return true
}

def checkContainers(containers, containersList) {
  containerExists = true
  containers.each{
    test = checkContainerExistence(it, containersList)
    !(test) ? containerExists = false : ""
  }
  return containerExists ? true : false
}

def checkParams(it) {
  // Check if params is in this given list
  return it in [
    'annotate-tools',
    'annotate-VCF',
    'annotateTools',
    'annotateVCF',
    'build',
    'call-name',
    'callName',
    'contact-mail',
    'contactMail',
    'container-path',
    'containerPath',
    'containers',
    'docker',
    'genome_base',
    'genome',
    'genomes',
    'help',
    'max_cpus',
    'max_memory',
    'max_time',
    'no-GVCF',
    'no-reports',
    'noGVCF',
    'noReports',
    'only-QC',
    'onlyQC',
    'out-dir',
    'outDir',
    'params',
    'project',
    'push',
    'repository',
    'sample-dir',
    'sample',
    'sampleDir',
    'single-CPUMem',
    'singleCPUMem',
    'singularity',
    'step',
    'tag',
    'test',
    'tools',
    'total-memory',
    'totalMemory',
    'vcflist',
    'verbose',
    'version']
}

def checkUppmaxProject() {
  // check if UPPMAX project number is specified
  return !(workflow.profile == 'slurm' && !params.project)
}

def defineContainersList(){
  // Return list of authorized containers
  return [
    'fastqc',
    'freebayes',
    'gatk',
    'igvtools',
    'multiqc',
    'mutect1',
    'picard',
    'qualimap',
    'r-base',
    'runallelecount',
    'sarek',
    'snpeff',
    'snpeffgrch37',
    'snpeffgrch38',
    'vepgrch37',
    'vepgrch38'
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
  log.info "       nextflow run SciLifeLab/Sarek/buildContainers.nf [--docker] [--push]"
  log.info "          [--containers <container1...>] [--singularity]"
  log.info "          [--containerPath <path>]"
  log.info "          [--tag <tag>] [--repository <repository>]"
  log.info "    Example:"
  log.info "      nextflow run SciLifeLab/Sarek/buildContainers.nf --docker --containers sarek"
  log.info "    --containers: Choose which containers to build"
  log.info "       Default: all"
  log.info "       Possible values:"
  log.info "         all, fastqc, freebayes, gatk, igvtools, multiqc, mutect1"
  log.info "         picard, qualimap, r-base, runallelecount, sarek"
  log.info "         snpeff, snpeffgrch37, snpeffgrch38, vepgrch37, vepgrch38"
  log.info "    --docker: Build containers using Docker"
  log.info "    --help"
  log.info "       you're reading it"
  log.info "    --push: Push containers to DockerHub"
  log.info "    --repository: Build containers under given repository"
  log.info "       Default: maxulysse"
  log.info "    --singularity: Download Singularity images"
  log.info "    --containerPath: Select where to download images"
  log.info "       Default: \$PWD"
  log.info "    --tag`: Choose the tag for the containers"
  log.info "       Default (version number): " + version
  log.info "    --version"
  log.info "       displays version number and more informations"
}

def isAllowedParams(params) {
  // Compare params to list of verified params
  final test = true
  params.each{
    if (!checkParams(it.toString().split('=')[0])) {
      println "params ${it.toString().split('=')[0]} is unknown"
      test = false
    }
  }
  return test
}

def minimalInformationMessage() {
  // Minimal information message
  log.info "Command Line: " + workflow.commandLine
  log.info "Project Dir : " + workflow.projectDir
  log.info "Launch Dir  : " + workflow.launchDir
  log.info "Work Dir    : " + workflow.workDir
  log.info "Containers  :"
  if (params.repository) log.info "  Repository   : ${params.repository}"
  else log.info "  ContainerPath: " + params.containerPath
  log.info "  Tag          : " + params.tag
}

def nextflowMessage() {
  // Nextflow message (version + build)
  log.info "N E X T F L O W  ~  version ${workflow.nextflow.version} ${workflow.nextflow.build}"
}

def startMessage() {
  // Display start message
  this.sarekMessage()
  this.minimalInformationMessage()
}

def versionMessage() {
  // Display version message
  log.info "CANCER ANALYSIS WORKFLOW"
  log.info "  version   : " + version
  log.info workflow.commitId ? "Git info    : ${workflow.repository} - ${workflow.revision} [${workflow.commitId}]" : "  revision  : " + this.grabRevision()
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
