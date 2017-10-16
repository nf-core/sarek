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
 - ProcessReference - Download all references if needed
 - DecompressFile - Extract files if needed
 - BuildBWAindexes - Build indexes for BWA
 - BuildPicardIndex - Build index with Picard
 - BuildSAMToolsIndex - Build index with SAMTools
 - BuildVCFIndex - Build index for VCF files
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

version = '1.2.3'

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
nf_required_version = '0.25.0'
try {
    if( ! nextflow.version.matches(">= $nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please update Nextflow.\n" +
              "============================================================"
}

if (params.help) exit 0, helpMessage()
if (params.version) exit 0, versionMessage()
if (!isAllowedParams(params)) exit 1, "params is unknown, see --help for more information"
if (!checkUppmaxProject()) exit 1, "No UPPMAX project ID found! Use --project <UPPMAX Project ID>"

// Default params:
// Such params are overridden by command line or configuration definitions

// No download
params.download = false
// refDir is empty
params.refDir = ''

verbose = params.verbose
download = params.download ? true : false

if (!download && params.refDir == "" ) exit 1, "No --refDir specified"
if (download && params.refDir != "" ) exit 1, "No need to specify --refDir"

if (params.genome == "smallGRCh37") {
  referencesFiles =
    [
      '1000G_phase1.indels.b37.small.vcf.gz',
      '1000G_phase3_20130502_SNP_maf0.3.small.loci',
      'b37_cosmic_v74.noCHR.sort.4.1.small.vcf.gz',
      'dbsnp_138.b37.small.vcf.gz',
      'human_g1k_v37_decoy.small.fasta.gz',
      'Mills_and_1000G_gold_standard.indels.b37.small.vcf.gz',
      'small.intervals'
    ]
} else if (params.genome == "GRCh37") {
  referencesFiles =
    [
      '1000G_phase1.indels.b37.vcf.gz',
      '1000G_phase3_20130502_SNP_maf0.3.loci.tar.bz2',
      'b37_cosmic_v74.noCHR.sort.4.1.vcf.tar.bz2',
      'dbsnp_138.b37.vcf.gz',
      'human_g1k_v37_decoy.fasta.gz',
      'Mills_and_1000G_gold_standard.indels.b37.vcf.gz',
      'wgs_calling_regions.grch37.list'
    ]
} else exit 1, "Can't build this reference genome"

if (download && params.genome != "smallGRCh37") exit 1, "Not possible to download $params.genome references files"

if (!download) referencesFiles.each{checkFile(params.refDir + "/" + it)}

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

startMessage()

process ProcessReference {
  tag download ? {"Download: " + reference} : {"Link: " + reference}

  input:
    val(reference) from referencesFiles

  output:
    file(reference) into processedFiles

  script:

  if (download)
  """
  wget https://github.com/szilvajuhos/smallRef/raw/master/$reference
  """

  else
  """
  ln -s $params.refDir/$reference .
  """
}


if (verbose) processedFiles = processedFiles.view {
  "Files preprocessed  : $it.fileName"
}

compressedfiles = Channel.create()
notCompressedfiles = Channel.create()

processedFiles
  .choice(compressedfiles, notCompressedfiles) {it =~ ".(gz|tar.bz2)" ? 0 : 1}

process DecompressFile {
  tag {reference}

  input:
    file(reference) from compressedfiles

  output:
    file("*.{vcf,fasta,loci}") into decompressedFiles

  script:
  realReference="readlink $reference"
  if (reference =~ ".gz")
    """
    gzip -d -c \$($realReference) > $reference.baseName
    """
  else if (reference =~ ".tar.bz2")
    """
    tar xvjf \$($realReference)
    """
}

if (verbose) decompressedFiles = decompressedFiles.view {
  "Files decomprecessed: $it.fileName"
}

fastaFile = Channel.create()
otherFiles = Channel.create()
vcfFiles = Channel.create()

decompressedFiles
  .choice(fastaFile, vcfFiles, otherFiles) {
    it =~ ".fasta" ? 0 :
    it =~ ".vcf" ? 1 : 2}

notCompressedfiles
  .mix(otherFiles)
  .collectFile(storeDir: "References/" + params.genome)

fastaForBWA = Channel.create()
fastaForPicard = Channel.create()
fastaForSAMTools = Channel.create()

fastaFile.into(fastaForBWA,fastaForPicard,fastaForSAMTools)

process BuildBWAindexes {
  tag {reference}

  publishDir "References/" + params.genome, mode: 'copy'

  input:
    file(reference) from fastaForBWA

  output:
    file(reference) into fastaFileToKeep
    file("*.{amb,ann,bwt,pac,sa}") into bwaIndexes

  script:

  """
  bwa index $reference
  """
}

if (verbose) fastaFileToKeep.view {
  "Fasta File          : $it.fileName"
}
if (verbose) bwaIndexes.flatten().view {
  "BWA index           : $it.fileName"
}

process BuildPicardIndex {
  tag {reference}

  publishDir "References/" + params.genome, mode: 'copy'

  input:
    file(reference) from fastaForPicard

  output:
    file("*.dict") into picardIndex

  script:
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar \$PICARD_HOME/picard.jar \
  CreateSequenceDictionary \
  REFERENCE=$reference \
  OUTPUT=${reference.baseName}.dict
  """
}

if (verbose) picardIndex.view {
  "Picard index        : $it.fileName"
}

process BuildSAMToolsIndex {
  tag {reference}

  publishDir "References/" + params.genome, mode: 'copy'

  input:
    file(reference) from fastaForSAMTools

  output:
    file("*.fai") into samtoolsIndex

  script:
  """
  samtools faidx $reference
  """
}

if (verbose) samtoolsIndex.view {
  "SAMTools index      : $it.fileName"
}

process BuildVCFIndex {
  tag {reference}

  publishDir "References/" + params.genome, mode: 'copy'

  input:
    file(reference) from vcfFiles

  output:
    file(reference) into vcfIndexed
    file("*.idx") into vcfIndex

  script:
  """
  \$IGVTOOLS_HOME/igvtools index $reference
  """
}

if (verbose) vcfIndexed.view {
  "VCF indexed         : $it.fileName"
}
if (verbose) vcfIndex.view {
  "VCF index           : $it.fileName"
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

def checkFile(it) {
  // Check file existence
  final f = file(it)
  if (!f.exists()) exit 1, "Missing file: $it, see --help for more information"
  return true
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
    'download',
    'genome',
    'genomes',
    'help',
    'no-GVCF',
    'no-reports',
    'noGVCF',
    'noReports',
    'project',
    'push',
    'ref-dir',
    'refDir',
    'repository',
    'run-time',
    'runTime',
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

def grabRevision() {
  // Return the same string executed from github or not
  return workflow.revision ?: workflow.commitId ?: workflow.scriptId.substring(0,10)
}

def helpMessage() {
  // Display help message
  this.cawMessage()
  log.info "    Usage:"
  log.info "       nextflow run buildReferences.nf --refDir <pathToRefDir> --genome <genome>"
  log.info "       nextflow run buildReferences.nf --download --genome smallGRCh37"
  log.info "       nextflow run SciLifeLab/CAW --test [--step STEP] [--tools TOOL[,TOOL]] --genome <Genome>"
  log.info "    --download"
  log.info "       Download reference files. (only with --genome smallGRCh37)"
  log.info "    --refDir <Directoy>"
  log.info "       Specify a directory containing reference files."
  log.info "    --genome <Genome>"
  log.info "       Use a specific genome version."
  log.info "       Possible values are:"
  log.info "         GRCh37"
  log.info "         smallGRCh37 (Build a small reference (for tests))"
  log.info "    --help"
  log.info "       you're reading it"
  log.info "    --version"
  log.info "       displays version number"
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
  log.info "Genome      : " + params.genome
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
