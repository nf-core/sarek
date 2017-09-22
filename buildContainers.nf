#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
================================================================================
=                          C A W - c o n ta i n e r s                          =
================================================================================
@Author
Maxime Garcia <maxime.garcia@scilifelab.se> [@MaxUlysse]
--------------------------------------------------------------------------------
 @Homepage
 https://github.com/SciLifeLab/CAW-containers
--------------------------------------------------------------------------------
 @Documentation
 https://github.com/SciLifeLab/CAW-containers/blob/master/README.md
--------------------------------------------------------------------------------
@Licence
 https://github.com/SciLifeLab/CAW-containers/blob/master/LICENSE
--------------------------------------------------------------------------------
 Processes overview
 - BuildContainers - Build containers using Docker
 - PushContainers - Push containers to DockerHub
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

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

if (params.help) {
  helpMessage()
  exit 1
}
if (params.version) {
  versionMessage()
  exit 1
}

startMessage()

if (!checkContainers(containers,containersList)) {exit 1, 'Unknown container(s), see --help for more information'}

/*
================================================================================
=                                 P R O C E S S                                =
================================================================================
*/

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

dockerContainersBuilt = dockerContainersBuilt.view {"Docker container: $repository/$it:$tag built."}

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

singularityContainersPulled = singularityContainersPulled.view {"Singularity container: $it pulled."}

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

dockerContainersPushed = dockerContainersPushed.view {"Docker container: $repository/$it:$tag pushed"}

/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/
def cawContainersMessage() {
  // Display CAW message
  log.info "CAW-containers ~ $version - " + this.grabRevision() + (workflow.commitId ? " [$workflow.commitId]" : "")
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

def defineContainersList(){
  // Return list of authorized containers
  return [
    'bcftools',
    'concatvcf',
    'fastqc',
    'freebayes',
    'gatk',
    'htslib',
    'igvtools',
    'mapreads',
    'multiqc',
    'mutect1',
    'picard',
    'qualimap',
    'runallelecount',
    'runascat',
    'runconvertallelecounts',
    'runmanta',
    'samtools',
    'snpeff',
    'snpeffgrch37',
    'snpeffgrch38',
    'strelka',
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
  this.cawContainersMessage()
  log.info "    Usage:"
  log.info "       nextflow run SciLifeLab/CAW-containers [--docker] [--push]"
  log.info "          [--containers <container1...>] [--singularity]"
  log.info "          [--singularityPublishDir <path>]"
  log.info "          [--tag <tag>] [--repository <repository>]"
  log.info "    Example:"
  log.info "      nextflow run . --docker --containers multiqc,fastqc"
  log.info "    --containers: Choose which containers to build"
  log.info "       Default: all"
  log.info "       Possible values:"
  log.info "         all, bcftools, concatvcf, fastqc, freebayes, gatk,"
  log.info "         htslib, igvtools, mapreads, multiqc, mutect1, picard,"
  log.info "         qualimap, runallelecount, runascat, runconvertallelecounts,"
  log.info "         runmanta, samtools, snpeff, snpeffgrch37, snpeffgrch38,"
  log.info "         strelka, vep, vepgrch37, vepgrch38"
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
  this.cawContainersMessage()
  this.minimalInformationMessage()
}

def versionMessage() {
  // Display version message
  log.info "CAW-containers"
  log.info "  version   : $version"
  log.info workflow.commitId ? "Git info    : $workflow.repository - $workflow.revision [$workflow.commitId]" : "  revision  : " + this.grabRevision()
}

workflow.onComplete {
  // Display complete message
  this.nextflowMessage()
  this.cawContainersMessage()
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
  this.cawContainersMessage()
  log.info "Workflow execution stopped with the following message: " + workflow.errorMessage
}
