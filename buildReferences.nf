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
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/


version = '1.1'

if (!isAllowedParams(params)) {exit 1, "params is unknown, see --help for more information"}

if (params.help) {
  helpMessage()
  exit 1
}

if (params.version) {
  versionMessage()
  exit 1
}

if (!checkUppmaxProject()) {exit 1, 'No UPPMAX project ID found! Use --project <UPPMAX Project ID>'}

params.download = false
params.refDir = ""

download = params.download ? true : false
if (!download && params.refDir == "" ) { exit 1, "No --refDir specified"}
if (download && params.refDir != "" ) { exit 1, "No need to specify --refDir"}

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
      'centromeres.list'
    ]
} else {exit 1, "Can't build this reference genome"}

if (download && params.genome != "smallGRCh37") {exit 1, "Not possible to download $params.genome references files"}

if (!download) {referencesFiles.each{checkFile(params.refDir + "/" + it)}}

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
  set -euo pipefail
  wget https://github.com/szilvajuhos/smallRef/raw/master/$reference
  """

  else
  """
  set -euo pipefail
  ln -s $params.refDir/$reference .
  """
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
   if (reference =~ ".gz")
     """
     set -euo pipefail
     realReference=`readlink $reference`
     gzip -d -c \$realReference > $reference.baseName
     """
   else if (reference =~ ".tar.bz2")
     """
     set -euo pipefail
     realReference=`readlink $reference`
     tar xvjf \$realReference
     """
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
    set file(reference), file("*.{amb,ann,bwt,pac,sa}") into bwaIndexes

  script:

  """
  set -euo pipefail
  bwa index $reference
  """
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
  set -euo pipefail
  java -Xmx${task.memory.toGiga()}g \
  -jar \$PICARD_HOME/picard.jar \
  CreateSequenceDictionary \
  REFERENCE=$reference \
  OUTPUT=${reference.baseName}.dict
  """
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
  set -euo pipefail
  samtools faidx $reference
  """
}

process BuildVCFIndex {
  tag {reference}

  publishDir "References/" + params.genome, mode: 'copy'

  input:
    file(reference) from vcfFiles

  output:
    set file(reference), file("*.idx") into vcfIndexed

  script:
  """
  set -euo pipefail
  \$IGVTOOLS_HOME/igvtools index $reference
  """
}

/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def checkFile(it) {
  // Check file existence
  final f = file(it)
  if (!f.exists()) {
    exit 1, "Missing file: $it, see --help for more information"
  }
  return true
}

def checkParams(it) {
  // Check if params is in this given list
  return it in [
    'annotate-tools',
    'annotate-VCF',
    'annotateTools',
    'annotateVCF',
    'call-name',
    'callName',
    'contact-mail',
    'contactMail',
    'download',
    'genome',
    'genomes',
    'help',
    'project',
    'ref-dir',
    'refDir',
    'reports',
    'run-time',
    'runTime',
    'sample-dir',
    'sample',
    'sampleDir',
    'single-CPUMem',
    'singleCPUMem',
    'step',
    'test',
    'tools',
    'vcflist',
    'verbose',
    'version']
}

def checkUppmaxProject() {
  return !((workflow.profile == 'standard' || workflow.profile == 'interactive') && !params.project)
}

def grabRevision() {
  return workflow.revision ?: workflow.scriptId.substring(0,10)
}

def helpMessage() { // Display help message
  log.info "CANCER ANALYSIS WORKFLOW ~ $version - revision: " + this.grabRevision()
  log.info "    Usage:"
  log.info "       nextflow run buildReferences.nf --refDir <pathToRefDir> --genome <genome>"
  log.info "       nextflow run buildReferences.nf --download --genome smallGRCh37"
  log.info "       nextflow run SciLifeLab/CAW --test [--step STEP] [--tools TOOL[,TOOL]] --genome <Genome>"
  log.info "    --download"
  log.info "       Download smallGRCh37 reference files."
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
  final test = true
  params.each{
    if (!checkParams(it.toString().split('=')[0])) {
      println "params ${it.toString().split('=')[0]} is unknown"
      test = false
    }
  }
  return test
}

def startMessage() { // Display start message
  log.info "CANCER ANALYSIS WORKFLOW ~ $version - revision: " + this.grabRevision()
  log.info "Command Line: $workflow.commandLine"
  log.info "Project Dir : $workflow.projectDir"
  log.info "Launch Dir  : $workflow.launchDir"
  log.info "Work Dir    : $workflow.workDir"
  log.info "Genome      : " + params.genome
}

def versionMessage() { // Display version message
  log.info "CANCER ANALYSIS WORKFLOW"
  log.info "  version   : $version"
  log.info workflow.commitId ? "Git info    : $workflow.repository - $workflow.revision [$workflow.commitId]" : "  revision  : " + this.grabRevision()
}

workflow.onComplete { // Display complete message
  log.info "N E X T F L O W ~ $workflow.nextflow.version - $workflow.nextflow.build"
  log.info "CANCER ANALYSIS WORKFLOW ~ $version - revision: " + this.grabRevision()
  log.info "Command Line: $workflow.commandLine"
  log.info "Project Dir : $workflow.projectDir"
  log.info "Launch Dir  : $workflow.launchDir"
  log.info "Work Dir    : $workflow.workDir"
  log.info "Genome      : " + params.genome
  log.info "Completed at: $workflow.complete"
  log.info "Duration    : $workflow.duration"
  log.info "Success     : $workflow.success"
  log.info "Exit status : $workflow.exitStatus"
  log.info "Error report: " + (workflow.errorReport ?: '-')
}

workflow.onError { // Display error message
  log.info "N E X T F L O W ~ version $workflow.nextflow.version [$workflow.nextflow.build]"
  log.info workflow.commitId ? "CANCER ANALYSIS WORKFLOW ~ $version - $workflow.revision [$workflow.commitId]" : "CANCER ANALYSIS WORKFLOW ~ $version - revision: " + this.grabRevision()
  log.info "Workflow execution stopped with the following message: " + workflow.errorMessage
}
