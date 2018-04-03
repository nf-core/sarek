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

if (!params.download && params.refDir == "" ) exit 1, "No --refDir specified"
if (params.download && params.refDir != "" ) exit 1, "No need to specify --refDir"

ch_referencesFiles = defReferencesFiles(params.genome)

if (params.download && params.genome != "smallGRCh37") exit 1, "Not possible to download ${params.genome} references files"

if (!params.download) ch_referencesFiles.each{checkFile(params.refDir + "/" + it)}

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

startMessage()

process ProcessReference {
  tag params.download ? {"Download: " + f_reference} : {"Link: " + f_reference}

  input:
    val(f_reference) from ch_referencesFiles

  output:
    file(f_reference) into ch_processedFiles

  script:

  if (params.download)
  """
  wget https://github.com/szilvajuhos/smallRef/raw/master/${f_reference}
  """

  else
  """
  ln -s ${params.refDir}/${f_reference} .
  """
}


if (params.verbose) ch_processedFiles = ch_processedFiles.view {
  "Files preprocessed  : ${it.fileName}"
}

ch_compressedfiles = Channel.create()
ch_notCompressedfiles = Channel.create()

ch_processedFiles
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
ch_otherFiles = Channel.create()
ch_vcfFiles = Channel.create()

ch_decompressedFiles
  .choice(ch_fastaFile, ch_vcfFiles, ch_otherFiles) {
    it =~ ".fasta" ? 0 :
    it =~ ".vcf" ? 1 : 2}

(ch_fastaFile, ch_fastaFileToKeep) = ch_fastaFile.into(2)
(ch_vcfFiles, ch_vcfFilesToKeep) = ch_vcfFiles.into(2)

ch_notCompressedfiles
  .mix(ch_otherFiles, ch_fastaFileToKeep, ch_vcfFilesToKeep)
  .collectFile(storeDir: params.outDir)

ch_fastaForBWA = Channel.create()
ch_fastaForPicard = Channel.create()
ch_fastaForSAMTools = Channel.create()

ch_fastaFile.into(ch_fastaForBWA,ch_fastaForPicard,ch_fastaForSAMTools)

process BuildBWAindexes {
  tag {f_reference}

  publishDir params.outDir, mode: 'link'

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

process BuildPicardIndex {
  tag {f_reference}

  publishDir params.outDir, mode: 'link'

  input:
    file(f_reference) from ch_fastaForPicard

  output:
    file("*.dict") into ch_picardIndex

  script:
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar \$PICARD_HOME/picard.jar \
  CreateSequenceDictionary \
  REFERENCE=${f_reference} \
  OUTPUT=${f_reference.baseName}.dict
  """
}

if (params.verbose) ch_picardIndex.view {
  "Picard index        : ${it.fileName}"
}

process BuildSAMToolsIndex {
  tag {f_reference}

  publishDir params.outDir, mode: 'link'

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

  publishDir params.outDir, mode: 'link'

  input:
    file(f_reference) from ch_vcfFiles

  output:
    file("${f_reference}.idx") into ch_vcfIndex

  script:
  """
  \$IGVTOOLS_HOME/igvtools index ${f_reference}
  """
}

if (params.verbose) ch_vcfIndex.view {
  "VCF index           : ${it.fileName}"
}

/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

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

def defReferencesFiles(genome) {
  if (genome == "smallGRCh37") {
    return [
    '1000G_phase1.indels.b37.small.vcf.gz',
    '1000G_phase3_20130502_SNP_maf0.3.small.loci',
    'b37_cosmic_v74.noCHR.sort.4.1.small.vcf.gz',
    'dbsnp_138.b37.small.vcf.gz',
    'human_g1k_v37_decoy.small.fasta.gz',
    'Mills_and_1000G_gold_standard.indels.b37.small.vcf.gz',
    'small.intervals'
    ]
  } else if (genome == "GRCh37") {
    return   [
    '1000G_phase1.indels.b37.vcf.gz',
    '1000G_phase3_20130502_SNP_maf0.3.loci.tar.bz2',
    'GRCh37_Cosmic_v83.vcf.tar.bz2',
    'dbsnp_138.b37.vcf.gz',
    'human_g1k_v37_decoy.fasta.gz',
    'Mills_and_1000G_gold_standard.indels.b37.vcf.gz',
    'wgs_calling_regions.grch37.list'
    ]
  } else exit 1, "Can't build this reference genome"
}

def grabRevision() {
  // Return the same string executed from github or not
  return workflow.revision ?: workflow.commitId ?: workflow.scriptId.substring(0,10)
}

def helpMessage() {
  // Display help message
  this.sarekMessage()
  log.info "    Usage:"
  log.info "       nextflow run buildReferences.nf --refDir <pathToRefDir> --genome <genome>"
  log.info "       nextflow run buildReferences.nf --download --genome smallGRCh37"
  log.info "    --download"
  log.info "       Download reference files. (only with --genome smallGRCh37)"
  log.info "    --refDir <Directoy>"
  log.info "       Specify a directory containing reference files."
  log.info "    --outDir <Directoy>"
  log.info "       Specify an output directory"
  log.info "    --genome <Genome>"
  log.info "       Choose which genome to build references from"
  log.info "       Possible values are:"
  log.info "         GRCh37"
  log.info "         smallGRCh37"
  log.info "    --help"
  log.info "       you're reading it"
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
  if (params.repository) log.info "  Repository   :" + params.repository
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
