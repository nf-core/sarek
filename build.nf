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
 Johannes Alneberg <johannes.alneberg@scilifelab.se> [@alneberg]
--------------------------------------------------------------------------------
 @Homepage
 https://sarek.scilifelab.se/
--------------------------------------------------------------------------------
 @Documentation
 https://github.com/SciLifeLab/Sarek/README.md
--------------------------------------------------------------------------------
 Processes overview
 - BuildWithDocker - Build containers using Docker
 - PullToSingularity - Pull Singularity containers from Docker Hub
 - PushToDocker - Push containers to Docker Hub
 - DecompressFile - Extract files if needed
 - BuildBWAindexes - Build indexes for BWA
 - BuildReferenceIndex - Build index for FASTA refs
 - BuildSAMToolsIndex - Build index with SAMTools
 - BuildVCFIndex - Build index for VCF files
 - BuildCache_snpEff - Download Cache for snpEff
 - BuildCache_VEP - Download taqbix index Cache for VEP
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
    if (!params.awsqueue) exit 1, "Provide the job queue for aws batch!"
}

// Define containers to handle (build/push or pull)
containersList = defineContainersList()
if (params.containers) {
  containers = params.containers.split(',').collect {it.trim()}
  containers = containers == ['all'] ? containersList : containers
} else containers = []

// push only to DockerHub, so only when using Docker
push = params.docker && params.push ? true : false

if (params.containers && !params.docker && !params.singularity) exit 1, 'No container technology choosed, specify --docker or --singularity, see --help for more information'

if (params.containers && !checkContainers(containers,containersList)) exit 1, 'Unknown container(s), see --help for more information'

ch_referencesFiles = Channel.fromPath("${params.refDir}/*")

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

startMessage()

/*
================================================================================
=                B  U  I  L  D      C  O  N  T  A  I  N  E  R  S               =
================================================================================
*/

dockerContainers = containers
singularityContainers = containers

process BuildWithDocker {
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

process PullToSingularity {
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

process PushToDocker {
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
=                B  U  I  L  D      R  E  F  E  R  E  N  C  E  S               =
================================================================================
*/

ch_compressedfiles = Channel.create()
ch_notCompressedfiles = Channel.create()

ch_referencesFiles
  .choice(ch_compressedfiles, ch_notCompressedfiles) {it =~ ".(gz|tar.bz2)" ? 0 : 1}

process DecompressFile {
  tag {f_reference}

  input:
    file(f_reference) from ch_compressedfiles

  output:
    file("*.{vcf,fasta,loci}") into ch_decompressedFiles

  script:
  realReferenceFile="readlink ${f_reference}"
  if (f_reference =~ ".gz")
    """
    gzip -d -c \$(${realReferenceFile}) > ${f_reference.baseName}
    """
  else if (f_reference =~ ".tar.bz2")
    """
    tar xvjf \$(${realReferenceFile})
    """
}

if (params.verbose) ch_decompressedFiles = ch_decompressedFiles.view {
  "Files decomprecessed: ${it.fileName}"
}

ch_fastaFile = Channel.create()
ch_fastaForBWA = Channel.create()
ch_fastaReference = Channel.create()
ch_fastaForSAMTools = Channel.create()
ch_otherFile = Channel.create()
ch_vcfFile = Channel.create()

ch_decompressedFiles
  .choice(ch_fastaFile, ch_vcfFile, ch_otherFile) {
    it =~ ".fasta" ? 0 :
    it =~ ".vcf" ? 1 : 2}

(ch_fastaForBWA, ch_fastaReference, ch_fastaForSAMTools, ch_fastaFileToKeep) = ch_fastaFile.into(4)
(ch_vcfFile, ch_vcfFileToKeep) = ch_vcfFile.into(2)

ch_notCompressedfiles
  .mix(ch_fastaFileToKeep, ch_vcfFileToKeep, ch_otherFile)
  .collectFile(storeDir: params.outDir)

process BuildBWAindexes {
  tag {f_reference}

  publishDir params.outDir, mode: params.publishDirMode

  input:
    file(f_reference) from ch_fastaForBWA

  output:
    file("*.{amb,ann,bwt,pac,sa}") into bwaIndexes

  script:
  """
  bwa index ${f_reference}
  """
}

if (params.verbose) bwaIndexes.flatten().view {
  "BWA index           : ${it.fileName}"
}

process BuildReferenceIndex {
  tag {f_reference}

  publishDir params.outDir, mode: params.publishDirMode

  input:
    file(f_reference) from ch_fastaReference

  output:
    file("*.dict") into ch_referenceIndex

  script:
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
  CreateSequenceDictionary \
  --REFERENCE ${f_reference} \
  --OUTPUT ${f_reference.baseName}.dict
  """
}

if (params.verbose) ch_referenceIndex.view {
  "Reference index     : ${it.fileName}"
}

process BuildSAMToolsIndex {
  tag {f_reference}

  publishDir params.outDir, mode: params.publishDirMode

  input:
    file(f_reference) from ch_fastaForSAMTools

  output:
    file("*.fai") into ch_samtoolsIndex

  script:
  """
  samtools faidx ${f_reference}
  """
}

if (params.verbose) ch_samtoolsIndex.view {
  "SAMTools index      : ${it.fileName}"
}

process BuildVCFIndex {
  tag {f_reference}

  publishDir params.outDir, mode: params.publishDirMode

  input:
    file(f_reference) from ch_vcfFile

  output:
    file("${f_reference}.idx") into ch_vcfIndex

  script:
  """
  igvtools index ${f_reference}
  """
}

if (params.verbose) ch_vcfIndex.view {
  "VCF index           : ${it.fileName}"
}

/*
================================================================================
=                   D  O  W  N  L  O  A  D      C  A  C  H  E                  =
================================================================================
*/

process BuildCache_snpEff {
  tag {snpeffDb}

  publishDir params.snpEff_cache, mode: params.publishDirMode

  input:
    val snpeffDb from Channel.value(params.genomes[params.genome].snpeffDb)

  output:
    file("*")

  when: params.snpEff_cache

  script:
  """
  snpEff download -v ${snpeffDb} -dataDir \${PWD}
  """
}

process BuildCache_VEP {
  tag {"${species}_${cache_version}_${genome}"}

  publishDir "${params.vep_cache}/${species}", mode: params.publishDirMode

  input:
    val cache_version from Channel.value(params.genomes[params.genome].vepCacheVersion)

  output:
    file("*")

  when: params.vep_cache

  script:
  genome = params.genome == "smallGRCh37" ? "GRCh37" : params.genome
  species = genome =~ "GRCh3*" ? "homo_sapiens" : ""
  """
  vep_install \
    -a cf \
    -c . \
    -s ${species} \
    -v ${cache_version} \
    -y ${genome} \
    --CACHE_VERSION ${cache_version} \
    --CONVERT \
    --NO_HTSLIB --NO_TEST --NO_BIOPERL --NO_UPDATE

  mv ${species}/* .
  rm -rf ${species}
  """
}

caddFileToDownload = (params.cadd_version) && (params.genome == "GRCh37" || params.genome == "GRCh38") ?
  Channel.from("https://krishna.gs.washington.edu/download/CADD/${params.cadd_version}/${params.genome}/InDels.tsv.gz",
    "https://krishna.gs.washington.edu/download/CADD/${params.cadd_version}/${params.genome}/whole_genome_SNVs.tsv.gz")
  : Channel.empty()

process DownloadCADD {
  tag {caddFile}

  publishDir "${params.cadd_cache}/${params.genome}", mode: params.publishDirMode

  input:
    val(caddFile) from caddFileToDownload

  output:
    set file("*.tsv.gz"), file("*.tsv.gz.tbi")

  when: params.cadd_cache

  script:
  """
  wget --quiet ${caddFile}
  tabix *.tsv.gz
  """
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

def checkFile(it) {
  // Check file existence
  final f = file(it)
  if (!f.exists()) exit 1, "Missing file: ${it}, see --help for more information"
  return true
}

def checkUppmaxProject() {
  // check if UPPMAX project number is specified
  return !(workflow.profile == 'slurm' && !params.project)
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

def grabRevision() {
  // Return the same string executed from github or not
  return workflow.revision ?: workflow.commitId ?: workflow.scriptId.substring(0,10)
}

def helpMessage() {
  // Display help message
  this.sarekMessage()
  log.info "    Usage:"
  log.info "    --help"
  log.info "       you're reading it"
  log.info "         BUILD CONTAINERS:"
  log.info "           nextflow run build.nf [--docker] [--push]"
  log.info "              [--containers <container1...>] [--singularity]"
  log.info "              [--containerPath <path>]"
  log.info "              [--tag <tag>] [--repository <repository>]"
  log.info "        --containers: Choose which containers to build"
  log.info "           Default: all"
  log.info "           Possible values:"
  log.info "             all, r-base, runallelecount, sarek"
  log.info "             snpeffgrch37, snpeffgrch38, vepgrch37, vepgrch38"
  log.info "        --docker: Build containers using Docker"
  log.info "        --push: Push containers to DockerHub"
  log.info "        --repository: Build containers under given repository"
  log.info "           Default: maxulysse"
  log.info "        --singularity: Download Singularity images"
  log.info "        --containerPath: Select where to download images"
  log.info "           Default: \$PWD"
  log.info "        --tag`: Choose the tag for the containers"
  log.info "           Default (version number): " + workflow.manifest.version
  log.info "         BUILD REFERENCES:"
  log.info "           nextflow run build.nf [--refDir <pathToRefDir> --outDir <pathToOutDir>]"
  log.info "        --refDir <Directoy>"
  log.info "           Specify a directory containing reference files"
  log.info "        --outDir <Directoy>"
  log.info "           Specify an output directory"
  log.info "         DOWNLOAD CACHE:"
  log.info "           nextflow run build.nf [--snpEff_cache <pathToSNPEFFcache>] [--vep_cache <pathToVEPcache>]"
  log.info "        --snpEff_cache <Directoy>"
  log.info "           Specify path to snpEff cache"
  log.info "           Will use snpEff version specified in configuration"
  log.info "        --vep_cache <Directoy>"
  log.info "           Specify path to VEP cache"
  log.info "           Will use VEP version specified in configuration"
  log.info "        --cadd_cache <Directoy>"
  log.info "           Specify path to CADD cache"
  log.info "           Will use CADD version specified"
  log.info "        --cadd_version <version>"
  log.info "           Will specify which CADD version to download"
}

def minimalInformationMessage() {
  // Minimal information message
  log.info "Command Line: " + workflow.commandLine
  log.info "Project Dir : " + workflow.projectDir
  log.info "Launch Dir  : " + workflow.launchDir
  log.info "Work Dir    : " + workflow.workDir
  log.info "Out Dir     : " + params.outDir
  log.info "Genome      : " + params.genome
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
