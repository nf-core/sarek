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
 - BuildDockerContainers - Build containers using Docker
 - PullSingularityContainers - Pull Singularity containers from Docker Hub
 - PushDockerContainers - Push containers to Docker Hub
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

// Define containers to handle (build/push or pull)
containersList = defineContainersList()
containers = params.containers.split(',').collect {it.trim()}
containers = containers == ['all'] ? containersList : containers

// push only to DockerHub, so only when using Docker
push = params.docker && params.push ? true : false

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
  tag {"${params.repository}/${container}:${params.tag}"}

  input:
    val container from dockerContainers

  output:
    val container into containersBuilt

  when: params.docker

  script:
  path = container == "sarek" ? "${baseDir}" : "${baseDir}/containers/${container}/."
  """
  docker build -t ${params.repository}/${container}:${params.tag} ${path}
  """
}

if (params.verbose) containersBuilt = containersBuilt.view {
  "Docker container: ${params.repository}/${it}:${params.tag} built."
}

process PullSingularityContainers {
  tag {"${params.repository}/${container}:${params.tag}"}

  publishDir "${params.containerPath}", mode: params.publishDirMode

  input:
    val container from singularityContainers

  output:
    file("${container}-${params.tag}.simg") into imagePulled

  when: params.singularity

  script:
  """
  singularity build ${container}-${params.tag}.simg docker://${params.repository}/${container}:${params.tag}
  """
}

if (params.verbose) imagePulled = imagePulled.view {
  "Singularity image: ${it.fileName} pulled."
}

process PushDockerContainers {
  tag {params.repository + "/" + container + ":" + params.tag}

  input:
    val container from containersBuilt

  output:
    val container into containersPushed

  when: params.docker && push

  script:
  """
  docker push ${params.repository}/${container}:${params.tag}
  """
}

if (params.verbose) containersPushed = containersPushed.view {
  "Docker container: ${params.repository}/${it}:${params.tag} pushed."
}

/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

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

def checkUppmaxProject() {
  // check if UPPMAX project number is specified
  return !(workflow.profile == 'slurm' && !params.project)
}

def grabRevision() {
  // Return the same string executed from github or not
  return workflow.revision ?: workflow.commitId ?: workflow.scriptId.substring(0,10)
}

def defineContainersList(){
  // Return list of authorized containers
  return [
    'r-base',
    'runallelecount',
    'sarek',
    'snpeffgrch37',
    'snpeffgrch38',
    'vepgrch37',
    'vepgrch38'
    ]
}

def helpMessage() {
  // Display help message
  this.sarekMessage()
  log.info "    Usage:"
  log.info "       nextflow run buildContainers.nf [--docker] [--push]"
  log.info "          [--containers <container1...>] [--singularity]"
  log.info "          [--containerPath <path>]"
  log.info "          [--tag <tag>] [--repository <repository>]"
  log.info "    Example:"
  log.info "      nextflow run buildContainers.nf --docker --containers sarek"
  log.info "    --containers: Choose which containers to build"
  log.info "       Default: all"
  log.info "       Possible values:"
  log.info "         all, r-base, runallelecount, sarek"
  log.info "         snpeffgrch37, snpeffgrch38, vepgrch37, vepgrch38"
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
  log.info "       Default (version number): " + workflow.manifest.version
}

def minimalInformationMessage() {
  // Minimal information message
  log.info "Command Line: " + workflow.commandLine
  log.info "Project Dir : " + workflow.projectDir
  log.info "Launch Dir  : " + workflow.launchDir
  log.info "Work Dir    : " + workflow.workDir
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
  log.info "Sarek ~ ${workflow.manifest.version} - " + this.grabRevision() + (workflow.commitId ? " [${workflow.commitId}]" : "")
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
