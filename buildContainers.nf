#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
kate: syntax groovy; space-indent on; indent-width 2;
================================================================================
=               C A N C E R    A N A L Y S I S    W O R K F L O W              =
================================================================================
 New Cancer Analysis Workflow. Started March 2016.
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
 http://opensource.scilifelab.se/projects/caw/
--------------------------------------------------------------------------------
 @Documentation
 https://github.com/SciLifeLab/CAW/README.md
--------------------------------------------------------------------------------
 Processes overview
 - BuildDockerContainers - Build containers using Docker
 - PullSingularityContainers - Pull Singularity containers from Docker Hub
 - PushDockerContainers - Push containers to Docker Hub
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

verbose = params.verbose
version = '1.1'
containersList = defineContainersList()
containers = params.containers.split(',').collect {it.trim()}
containers = containers == ['all'] ? containersList : containers
docker = params.docker ? true : false
push = params.docker && params.push ? true : false
repository = params.repository
tag = params.tag ? params.tag : version
singularity = params.singularity ? true : false
singularityPublishDir = params.singularity && params.singularityPublishDir ? params.singularityPublishDir : "."

if (params.help) exit 1, helpMessage()
if (params.version) exit 1, versionMessage()
if (!isAllowedParams(params)) exit 1, "params is unknown, see --help for more information"
if (!nextflow.version.matches('>= 0.25.0')) exit 1, "Nextflow version 0.25.0 or greater is needed to run this workflow"

if (!checkUppmaxProject()) exit 1, "No UPPMAX project ID found! Use --project <UPPMAX Project ID>"

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
  tag {repository + "/" + container + ":" + tag}

  input:
    val container from dockerContainers

  output:
    val container into dockerContainersBuilt

  when: docker

  script:
  """
  docker build -t $repository/$container:$tag $baseDir/containers/$container/.
  """
}

if (verbose) dockerContainersBuilt = dockerContainersBuilt.view {
  "Docker container: $repository/$it:$tag built."
}

process PullSingularityContainers {
  tag {repository + "/" + container + ":" + tag}

  publishDir singularityPublishDir, mode: 'move'

  input:
    val container from singularityContainers

  output:
    file("*.img") into singularityContainersPulled

  when: singularity

  script:
  """
  singularity pull --name $container-${tag}.img docker://$repository/$container:$tag
  """
}

if (verbose) singularityContainersPulled = singularityContainersPulled.view {
  "Singularity container: $it pulled."
}

process PushDockerContainers {
  tag {repository + "/" + container + ":" + tag}

  input:
    val container from dockerContainersBuilt

  output:
    val container into dockerContainersPushed

  when: docker && push

  script:
  """
  docker push $repository/$container:$tag
  """
}

if (verbose) dockerContainersPushed = dockerContainersPushed.view {
  "Docker container: $repository/$it:$tag pushed."
}

/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def cawMessage() {
  // Display CAW message
  log.info "CANCER ANALYSIS WORKFLOW ~ $version - " + this.grabRevision() + (workflow.commitId ? " [$workflow.commitId]" : "")
}

def checkContainerExistence(container, list) {
  try {assert list.contains(container)}
  catch (AssertionError ae) {
    println("Unknown container: $container")
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
    'containers',
    'docker',
    'genome',
    'genomes',
    'help',
    'no-GVCF',
    'no-reports',
    'noGVCF',
    'noReports',
    'project',
    'push',
    'repository',
    'sample-dir',
    'sample',
    'sampleDir',
    'single-CPUMem',
    'singleCPUMem',
    'singularity-publish-dir',
    'singularity',
    'singularityPublishDir',
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
    'caw',
    'fastqc',
    'freebayes',
    'gatk',
    'igvtools',
    'multiqc',
    'mutect1',
    'picard',
    'qualimap',
    'runallelecount',
    'runascat',
    'runconvertallelecounts',
    'snpeff',
    'snpeffgrch37',
    'snpeffgrch38',
    'vep',
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
  this.cawMessage()
  log.info "    Usage:"
  log.info "       nextflow run SciLifeLab/buildContainers.nf [--docker] [--push]"
  log.info "          [--containers <container1...>] [--singularity]"
  log.info "          [--singularityPublishDir <path>]"
  log.info "          [--tag <tag>] [--repository <repository>]"
  log.info "    Example:"
  log.info "      nextflow run . --docker --containers multiqc,fastqc"
  log.info "    --containers: Choose which containers to build"
  log.info "       Default: all"
  log.info "       Possible values:"
  log.info "         all, caw, fastqc, freebayes, gatk, igvtools, multiqc"
  log.info "         mutect1, picard, qualimap, runallelecount, runascat"
  log.info "         runconvertallelecounts, snpeff, snpeffgrch37, snpeffgrch38"
  log.info "         vep, vepgrch37, vepgrch38"
  log.info "    --docker: Build containers using Docker"
  log.info "    --help"
  log.info "       you're reading it"
  log.info "    --push: Push containers to DockerHub"
  log.info "    --repository: Build containers under given repository"
  log.info "       Default: maxulysse"
  log.info "    --singularity: Build containers using Singularity"
  log.info "    --singularityPublishDir: Select where to download containers"
  log.info "       Default: $PWD"
  log.info "    --tag`: Build containers using given tag"
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
  log.info "Command Line: $workflow.commandLine"
  log.info "Project Dir : $workflow.projectDir"
  log.info "Launch Dir  : $workflow.launchDir"
  log.info "Work Dir    : $workflow.workDir"
  log.info "Containers  : " + containers.join(', ')
}

def nextflowMessage() {
  // Nextflow message (version + build)
  log.info "N E X T F L O W  ~  version $workflow.nextflow.version $workflow.nextflow.build"
}

def startMessage() {
  // Display start message
  this.cawMessage()
  this.minimalInformationMessage()
}

def versionMessage() {
  // Display version message
  log.info "CANCER ANALYSIS WORKFLOW"
  log.info "  version   : $version"
  log.info workflow.commitId ? "Git info    : $workflow.repository - $workflow.revision [$workflow.commitId]" : "  revision  : " + this.grabRevision()
}

workflow.onComplete {
  // Display complete message
  this.nextflowMessage()
  this.cawMessage()
  this.minimalInformationMessage()
  log.info "Completed at: $workflow.complete"
  log.info "Duration    : $workflow.duration"
  log.info "Success     : $workflow.success"
  log.info "Exit status : $workflow.exitStatus"
  log.info "Error report: " + (workflow.errorReport ?: '-')
}

workflow.onError {
  // Display error message
  this.nextflowMessage()
  this.cawMessage()
  log.info "Workflow execution stopped with the following message: " + workflow.errorMessage
}
