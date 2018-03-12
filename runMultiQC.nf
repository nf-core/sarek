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
 - RunMultiQC - Run MultiQC on reports
 - MapReads - Map reads with BWA
 - MergeBams - Merge BAMs if multilane samples
 - MarkDuplicates - Mark Duplicates with Picard
 - RealignerTargetCreator - Create realignment target intervals
 - IndelRealigner - Realign BAMs as T/N pair
 - CreateRecalibrationTable - Create Recalibration Table with BaseRecalibrator
 - RecalibrateBam - Recalibrate Bam with PrintReads
 - RunSamtoolsStats - Run Samtools stats on recalibrated BAM files
 - RunBamQC - Run qualimap BamQC on recalibrated BAM files
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

version = '1.3'

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

// Reports are generated
params.noReports = false
// outDir is current directory
params.outDir = baseDir
// Params are defined in config files
params.containerPath = ''
params.repository = ''
params.tag = ''

directoryMap = defineDirectoryMap()

reports = !params.noReports
verbose = params.verbose

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

startMessage()

process GenerateMultiQCconfig {
  publishDir "${params.outDir}/${directoryMap.multiQC}", mode: 'copy'

  input:

  output:
  file("multiqc_config.yaml") into multiQCconfig

  when: reports

  script:
  """
  touch multiqc_config.yaml
  echo "custom_logo: ${baseDir}/doc/images/Sarek_no_Border.png" >> multiqc_config.yaml
  echo "custom_logo_url: http://opensource.scilifelab.se/projects/sarek" >> multiqc_config.yaml
  echo "custom_logo_title: 'Sarek'" >> multiqc_config.yaml
  echo "report_header_info:" >> multiqc_config.yaml
  echo "- Sarek version: ${version}" >> multiqc_config.yaml
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

if (verbose && reports) multiQCconfig = multiQCconfig.view {
  "MultiQC config:\n\
  File  : [${it.fileName}]"
}

reportsForMultiQC = Channel.empty()
  .mix(
    Channel.fromPath("${params.outDir}/Reports/bamQC/*", type: 'dir'),
    Channel.fromPath("${params.outDir}/Reports/BCFToolsStats/*"),
    Channel.fromPath("${params.outDir}/Reports/FastQC/*/*"),
    Channel.fromPath("${params.outDir}/Reports/MarkDuplicates/*"),
    Channel.fromPath("${params.outDir}/Reports/SamToolsStats/*"),
    multiQCconfig
  ).collect()

process RunMultiQC {
  publishDir "${params.outDir}/${directoryMap.multiQC}", mode: 'copy'

  input:
    file ('*') from reportsForMultiQC

  output:
    set file("*multiqc_report.html"), file("*multiqc_data") into multiQCReport

  when: reports

  script:
  """
  multiqc -f -v .
  """
}

if (verbose) multiQCReport = multiQCReport.view {
  "MultiQC report:\n\
  File  : [${it[0].fileName}]\n\
  Dir   : [${it[1].fileName}]"
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

def checkParameterExistence(it, list) {
  // Check parameter existence
  if (!list.contains(it)) {
    println("Unknown parameter: ${it}")
    return false
  }
  return true
}

def checkParams(it) {
  // Check if params is in this given list
  return it in [
    'ac-loci',
    'acLoci',
    'annotate-tools',
    'annotate-VCF',
    'annotateTools',
    'annotateVCF',
    'build',
    'bwa-index',
    'bwaIndex',
    'call-name',
    'callName',
    'contact-mail',
    'contactMail',
    'container-path',
    'containerPath',
    'containers',
    'cosmic-index',
    'cosmic',
    'cosmicIndex',
    'dbsnp-index',
    'dbsnp',
    'docker',
    'explicitBqsrNeeded',
    'genome_base',
    'genome-dict',
    'genome-file',
    'genome-index',
    'genome',
    'genomeDict',
    'genomeFile',
    'genomeIndex',
    'genomes',
    'help',
    'intervals',
    'known-indels-index',
    'known-indels',
    'knownIndels',
    'knownIndelsIndex',
    'max_cpus',
    'max_memory',
    'max_time',
    'no-BAMQC',
    'no-GVCF',
    'no-reports',
    'noBAMQC',
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

def defineDirectoryMap() {
  return [
    'multiQC'          : 'Reports/MultiQC'
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
  log.info "       nextflow run SciLifeLab/Sarek --sample <file.tsv> [--step STEP] --genome <Genome>"
  log.info "       nextflow run SciLifeLab/Sarek --sampleDir <Directory> [--step STEP] --genome <Genome>"
  log.info "       nextflow run SciLifeLab/Sarek --test [--step STEP] --genome <Genome>"
  log.info "    --sample <file.tsv>"
  log.info "       Specify a TSV file containing paths to sample files."
  log.info "    --sampleDir <Directoy>"
  log.info "       Specify a directory containing sample files."
  log.info "    --test"
  log.info "       Use a test sample."
  log.info "    --step"
  log.info "       Option to start workflow"
  log.info "       Possible values are:"
  log.info "         mapping (default, will start workflow with FASTQ files)"
  log.info "         realign (will start workflow with non-realigned BAM files)"
  log.info "         recalibrate (will start workflow with non-recalibrated BAM files)"
  log.info "    --noReports"
  log.info "       Disable QC tools and MultiQC to generate a HTML report"
  log.info "    --genome <Genome>"
  log.info "       Use a specific genome version."
  log.info "       Possible values are:"
  log.info "         GRCh37"
  log.info "         GRCh38 (Default)"
  log.info "         smallGRCh37 (Use a small reference (Tests only))"
  log.info "    --onlyQC"
  log.info "       Run only QC tools and gather reports"
  log.info "    --help"
  log.info "       you're reading it"
  log.info "    --verbose"
  log.info "       Adds more verbosity to workflow"
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

def startMessage() {
  // Display start message
  this.sarekMessage()
  this.minimalInformationMessage()
}

def versionMessage() {
  // Display version message
  log.info "Sarek"
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
