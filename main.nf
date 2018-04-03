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
 - RunFastQC - Run FastQC for QC on fastq files
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

step = params.step.toLowerCase()
if (step == 'preprocessing') step = 'mapping'

directoryMap = defineDirectoryMap()
referenceMap = defineReferenceMap()
stepList = defineStepList()

if (!checkParameterExistence(step, stepList)) exit 1, 'Unknown step, see --help for more information'
if (step.contains(',')) exit 1, 'You can choose only one step, see --help for more information'
if (step == 'mapping' && !checkExactlyOne([params.test, params.sample, params.sampleDir]))
  exit 1, 'Please define which samples to work on by providing exactly one of the --test, --sample or --sampleDir options'
if (!checkReferenceMap(referenceMap)) exit 1, 'Missing Reference file(s), see --help for more information'

if (params.test && params.genome in ['GRCh37', 'GRCh38']) {
  referenceMap.intervals = file("$workflow.projectDir/repeats/tiny_${params.genome}.list")
}

tsvPath = ''
if (params.sample) tsvPath = params.sample

// No need for tsv file for step annotate
if (!params.sample && !params.sampleDir) {
  tsvPaths = [
      'mapping': "${workflow.projectDir}/data/tsv/tiny.tsv",
      'realign': "${directoryMap.nonRealigned}/nonRealigned.tsv",
      'recalibrate': "${directoryMap.nonRecalibrated}/nonRecalibrated.tsv"
  ]
  if (params.test || step != 'mapping') tsvPath = tsvPaths[step]
}

// Set up the fastqFiles and bamFiles channels. One of them remains empty
// Except for step annotate, in which both stay empty
fastqFiles = Channel.empty()
bamFiles = Channel.empty()
if (tsvPath) {
  tsvFile = file(tsvPath)
  switch (step) {
    case 'mapping': fastqFiles = extractFastq(tsvFile); break
    case 'realign': bamFiles = extractBams(tsvFile); break
    case 'recalibrate': bamFiles = extractRecal(tsvFile); break
    default: exit 1, "Unknown step ${step}"
  }
} else if (params.sampleDir) {
  if (step != 'mapping') exit 1, '--sampleDir does not support steps other than "mapping"'
  fastqFiles = extractFastqFromDir(params.sampleDir)
  (fastqFiles, fastqTmp) = fastqFiles.into(2)
  fastqTmp.toList().subscribe onNext: {
    if (it.size() == 0) {
      exit 1, "No FASTQ files found in --sampleDir directory '${params.sampleDir}'"
    }
  }
  tsvFile = params.sampleDir  // used in the reports
} else exit 1, 'No sample were defined, see --help'

if (step == 'mapping') (patientGenders, fastqFiles) = extractGenders(fastqFiles)
else (patientGenders, bamFiles) = extractGenders(bamFiles)

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

startMessage()

(fastqFiles, fastqFilesforFastQC) = fastqFiles.into(2)

if (params.verbose) fastqFiles = fastqFiles.view {
  "FASTQs to preprocess:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\tRun   : ${it[3]}\n\
  Files : [${it[4].fileName}, ${it[5].fileName}]"
}

if (params.verbose) bamFiles = bamFiles.view {
  "BAMs to process:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  Files : [${it[3].fileName}, ${it[4].fileName}]"
}

process RunFastQC {
  tag {idPatient + "-" + idRun}

  publishDir "${directoryMap.fastQC}/${idRun}", mode: 'link'

  input:
    set idPatient, status, idSample, idRun, file(fastqFile1), file(fastqFile2) from fastqFilesforFastQC

  output:
    file "*_fastqc.{zip,html}" into fastQCreport

  when: step == 'mapping' && !params.noReports

  script:
  """
  fastqc -q ${fastqFile1} ${fastqFile2}
  """
}

if (params.verbose) fastQCreport = fastQCreport.view {
  "FastQC report:\n\
  Files : [${it[0].fileName}, ${it[1].fileName}]"
}

process MapReads {
  tag {idPatient + "-" + idRun}

  input:
    set idPatient, status, idSample, idRun, file(fastqFile1), file(fastqFile2) from fastqFiles
    set file(genomeFile), file(bwaIndex) from Channel.value([referenceMap.genomeFile, referenceMap.bwaIndex])

  output:
    set idPatient, status, idSample, idRun, file("${idRun}.bam") into mappedBam

  when: step == 'mapping' && !params.onlyQC

  script:
  readGroup = "@RG\\tID:${idRun}\\tPU:${idRun}\\tSM:${idSample}\\tLB:${idSample}\\tPL:illumina"
  // adjust mismatch penalty for tumor samples
  extra = status == 1 ? "-B 3 " : ""
  """
  bwa mem -R \"${readGroup}\" ${extra}-t ${task.cpus} -M \
  ${genomeFile} ${fastqFile1} ${fastqFile2} | \
  samtools sort --threads ${task.cpus} -m 4G - > ${idRun}.bam
  """
}

if (params.verbose) mappedBam = mappedBam.view {
  "Mapped BAM (single or to be merged):\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\tRun   : ${it[3]}\n\
  File  : [${it[4].fileName}]"
}

// Sort bam whether they are standalone or should be merged
// Borrowed code from https://github.com/guigolab/chip-nf

singleBam = Channel.create()
groupedBam = Channel.create()
mappedBam.groupTuple(by:[0,1,2])
  .choice(singleBam, groupedBam) {it[3].size() > 1 ? 1 : 0}
singleBam = singleBam.map {
  idPatient, status, idSample, idRun, bam ->
  [idPatient, status, idSample, bam]
}

process MergeBams {
  tag {idPatient + "-" + idSample}

  input:
    set idPatient, status, idSample, idRun, file(bam) from groupedBam

  output:
    set idPatient, status, idSample, file("${idSample}.bam") into mergedBam

  when: step == 'mapping' && !params.onlyQC

  script:
  """
  samtools merge --threads ${task.cpus} ${idSample}.bam ${bam}
  """
}

if (params.verbose) singleBam = singleBam.view {
  "Single BAM:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  File  : [${it[3].fileName}]"
}

if (params.verbose) mergedBam = mergedBam.view {
  "Merged BAM:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  File  : [${it[3].fileName}]"
}

mergedBam = mergedBam.mix(singleBam)

if (params.verbose) mergedBam = mergedBam.view {
  "BAM for MarkDuplicates:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  File  : [${it[3].fileName}]"
}

process MarkDuplicates {
  tag {idPatient + "-" + idSample}

  publishDir params.outDir, saveAs: { it == "${bam}.metrics" ? "${directoryMap.markDuplicatesQC}/${it}" : "${directoryMap.nonRealigned}/${it}" }, mode: 'link'

  input:
    set idPatient, status, idSample, file(bam) from mergedBam

  output:
    set idPatient, file("${idSample}_${status}.md.bam"), file("${idSample}_${status}.md.bai") into duplicates
    set idPatient, status, idSample, val("${idSample}_${status}.md.bam"), val("${idSample}_${status}.md.bai") into markDuplicatesTSV
    file ("${bam}.metrics") into markDuplicatesReport

  when: step == 'mapping' && !params.onlyQC

  script:
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar \$PICARD_HOME/picard.jar MarkDuplicates \
  INPUT=${bam} \
  METRICS_FILE=${bam}.metrics \
  TMP_DIR=. \
  ASSUME_SORTED=true \
  VALIDATION_STRINGENCY=LENIENT \
  CREATE_INDEX=TRUE \
  OUTPUT=${idSample}_${status}.md.bam
  """
}

// Creating a TSV file to restart from this step
markDuplicatesTSV.map { idPatient, status, idSample, bam, bai ->
  gender = patientGenders[idPatient]
  "${idPatient}\t${gender}\t${status}\t${idSample}\t${directoryMap.nonRealigned}/${bam}\t${directoryMap.nonRealigned}/${bai}\n"
}.collectFile(
  name: 'nonRealigned.tsv', sort: true, storeDir: "${directoryMap.nonRealigned}"
)

// Create intervals for realignement using both tumor+normal as input
// Group the marked duplicates BAMs for intervals and realign by idPatient
// Grouping also by gender, to make a nicer channel
duplicatesGrouped = Channel.empty()
if (step == 'mapping') duplicatesGrouped = duplicates.groupTuple()
else if (step == 'realign') duplicatesGrouped = bamFiles.map{
  idPatient, status, idSample, bam, bai ->
  [idPatient, bam, bai]
}.groupTuple()

// The duplicatesGrouped channel is duplicated
// one copy goes to the RealignerTargetCreator process
// and the other to the IndelRealigner process
(duplicatesInterval, duplicatesRealign) = duplicatesGrouped.into(2)

if (params.verbose) duplicatesInterval = duplicatesInterval.view {
  "BAMs for RealignerTargetCreator:\n\
  ID    : ${it[0]}\n\
  Files : ${it[1].fileName}\n\
  Files : ${it[2].fileName}"
}

if (params.verbose) duplicatesRealign = duplicatesRealign.view {
  "BAMs to phase:\n\
  ID    : ${it[0]}\n\
  Files : ${it[1].fileName}\n\
  Files : ${it[2].fileName}"
}

if (params.verbose) markDuplicatesReport = markDuplicatesReport.view {
  "MarkDuplicates report:\n\
  File  : [${it.fileName}]"
}

// VCF indexes are added so they will be linked, and not re-created on the fly
//  -L "1:131941-141339" \

process RealignerTargetCreator {
  tag {idPatient}

  input:
    set idPatient, file(bam), file(bai) from duplicatesInterval
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(knownIndels), file(knownIndelsIndex), file(intervals) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.knownIndels,
      referenceMap.knownIndelsIndex,
      referenceMap.intervals
    ])

  output:
    set idPatient, file("${idPatient}.intervals") into intervals

  when: ( step == 'mapping' || step == 'realign' ) && !params.onlyQC

  script:
  bams = bam.collect{"-I ${it}"}.join(' ')
  known = knownIndels.collect{"-known ${it}"}.join(' ')
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar \$GATK_HOME/GenomeAnalysisTK.jar \
  -T RealignerTargetCreator \
  ${bams} \
  -R ${genomeFile} \
  ${known} \
  -nt ${task.cpus} \
  -L ${intervals} \
  -o ${idPatient}.intervals
  """
}

if (params.verbose) intervals = intervals.view {
  "Intervals to phase:\n\
  ID    : ${it[0]}\n\
  File  : [${it[1].fileName}]"
}

bamsAndIntervals = duplicatesRealign
  .phase(intervals)
  .map{duplicatesRealign, intervals ->
    tuple(
      duplicatesRealign[0],
      duplicatesRealign[1],
      duplicatesRealign[2],
      intervals[1]
    )}

if (params.verbose) bamsAndIntervals = bamsAndIntervals.view {
  "BAMs and Intervals phased for IndelRealigner:\n\
  ID    : ${it[0]}\n\
  Files : ${it[1].fileName}\n\
  Files : ${it[2].fileName}\n\
  File  : [${it[3].fileName}]"
}

// use nWayOut to split into T/N pair again
process IndelRealigner {
  tag {idPatient}

  publishDir directoryMap.nonRecalibrated, mode: 'link'

  input:
    set idPatient, file(bam), file(bai), file(intervals) from bamsAndIntervals
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(knownIndels), file(knownIndelsIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.knownIndels,
      referenceMap.knownIndelsIndex])

  output:
    set idPatient, file("*.real.bam"), file("*.real.bai") into realignedBam mode flatten

  when: ( step == 'mapping' || step == 'realign' ) && !params.onlyQC

  script:
  bams = bam.collect{"-I ${it}"}.join(' ')
  known = knownIndels.collect{"-known ${it}"}.join(' ')
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar \$GATK_HOME/GenomeAnalysisTK.jar \
  -T IndelRealigner \
  ${bams} \
  -R ${genomeFile} \
  -targetIntervals ${intervals} \
  ${known} \
  -nWayOut '.real.bam'
  """
}

realignedBam = realignedBam.map {
    idPatient, bam, bai ->
    tag = bam.baseName.tokenize('.')[0]
    status   = tag[-1..-1].toInteger()
    idSample = tag.take(tag.length()-2)
    [idPatient, status, idSample, bam, bai]
}

if (params.verbose) realignedBam = realignedBam.view {
  "Realigned BAM to CreateRecalibrationTable:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  Files : [${it[3].fileName}, ${it[4].fileName}]"
}

process CreateRecalibrationTable {
  tag {idPatient + "-" + idSample}

  publishDir directoryMap.nonRecalibrated, mode: 'link', overwrite: false

  input:
    set idPatient, status, idSample, file(bam), file(bai) from realignedBam
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(dbsnp), file(dbsnpIndex), file(knownIndels), file(knownIndelsIndex), file(intervals) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.dbsnp,
      referenceMap.dbsnpIndex,
      referenceMap.knownIndels,
      referenceMap.knownIndelsIndex,
      referenceMap.intervals,
    ])

  output:
    set idPatient, status, idSample, file(bam), file(bai), file("${idSample}.recal.table") into recalibrationTable
    set idPatient, status, idSample, val("${idSample}_${status}.md.real.bam"), val("${idSample}_${status}.md.real.bai"), val("${idSample}.recal.table") into recalibrationTableTSV

  when: ( step == 'mapping' || step == 'realign' ) && !params.onlyQC

  script:
  known = knownIndels.collect{ "-knownSites ${it}" }.join(' ')
  """
  java -Xmx${task.memory.toGiga()}g \
  -Djava.io.tmpdir="/tmp" \
  -jar \$GATK_HOME/GenomeAnalysisTK.jar \
  -T BaseRecalibrator \
  -R ${genomeFile} \
  -I ${bam} \
  -L ${intervals} \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  -knownSites ${dbsnp} \
  ${known} \
  -nct ${task.cpus} \
  -l INFO \
  -o ${idSample}.recal.table
  """
}
// Create a TSV file to restart from this step
recalibrationTableTSV.map { idPatient, status, idSample, bam, bai, recalTable ->
  gender = patientGenders[idPatient]
  "${idPatient}\t${gender}\t${status}\t${idSample}\t${directoryMap.nonRecalibrated}/${bam}\t${directoryMap.nonRecalibrated}/${bai}\t${directoryMap.nonRecalibrated}/${recalTable}\n"
}.collectFile(
  name: 'nonRecalibrated.tsv', sort: true, storeDir: directoryMap.nonRecalibrated
)

if (step == 'recalibrate') recalibrationTable = bamFiles

if (params.verbose) recalibrationTable = recalibrationTable.view {
  "Base recalibrated table for RecalibrateBam:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  Files : [${it[3].fileName}, ${it[4].fileName}, ${it[5].fileName}]"
}

(bamForBamQC, bamForSamToolsStats, recalTables, recalibrationTableForHC, recalibrationTable) = recalibrationTable.into(5)

// Remove recalTable from Channels to match inputs for Process to avoid:
// WARN: Input tuple does not match input set cardinality declared by process...
bamForBamQC = bamForBamQC.map { it[0..4] }
bamForSamToolsStats = bamForSamToolsStats.map{ it[0..4] }

recalTables = recalTables.map { [it[0]] + it[2..-1] } // remove status

process RecalibrateBam {
  tag {idPatient + "-" + idSample}

  publishDir directoryMap.recalibrated, mode: 'link'

  input:
    set idPatient, status, idSample, file(bam), file(bai), file(recalibrationReport) from recalibrationTable
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(intervals) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.intervals,
    ])

  output:
    set idPatient, status, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bai") into recalibratedBam, recalibratedBamForStats
    set idPatient, status, idSample, val("${idSample}.recal.bam"), val("${idSample}.recal.bai") into recalibratedBamTSV

  // HaplotypeCaller can do BQSR on the fly, so do not create a
  // recalibrated BAM explicitly.
  when: params.explicitBqsrNeeded && !params.onlyQC

  script:
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar \$GATK_HOME/GenomeAnalysisTK.jar \
  -T PrintReads \
  -R ${genomeFile} \
  -I ${bam} \
  -L ${intervals} \
  --BQSR ${recalibrationReport} \
  -o ${idSample}.recal.bam
  """
}
// Creating a TSV file to restart from this step
recalibratedBamTSV.map { idPatient, status, idSample, bam, bai ->
  gender = patientGenders[idPatient]
  "${idPatient}\t${gender}\t${status}\t${idSample}\t${directoryMap.recalibrated}/${bam}\t${directoryMap.recalibrated}/${bai}\n"
}.collectFile(
  name: 'recalibrated.tsv', sort: true, storeDir: directoryMap.recalibrated
)

if (params.verbose) recalibratedBam = recalibratedBam.view {
  "Recalibrated BAM for variant Calling:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  Files : [${it[3].fileName}, ${it[4].fileName}]"
}

process RunSamtoolsStats {
  tag {idPatient + "-" + idSample}

  publishDir directoryMap.samtoolsStats, mode: 'link'

  input:
    set idPatient, status, idSample, file(bam), file(bai) from bamForSamToolsStats

  output:
    file ("${bam}.samtools.stats.out") into samtoolsStatsReport

  when: !params.noReports

  script:
  """
  samtools stats ${bam} > ${bam}.samtools.stats.out
  """
}

if (params.verbose) samtoolsStatsReport = samtoolsStatsReport.view {
  "SAMTools stats report:\n\
  File  : [${it.fileName}]"
}

process RunBamQC {
  tag {idPatient + "-" + idSample}

  publishDir directoryMap.bamQC, mode: 'link'

  input:
    set idPatient, status, idSample, file(bam), file(bai) from bamForBamQC

  output:
    file("${idSample}") into bamQCreport

  when: !params.noReports && !params.noBAMQC

  script:
  """
  qualimap --java-mem-size=${task.memory.toGiga()}G \
  bamqc \
  -bam ${bam} \
  -outdir ${idSample} \
  -outformat HTML
  """
}

if (params.verbose) bamQCreport = bamQCreport.view {
  "BamQC report:\n\
  Dir   : [${it.fileName}]"
}

/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def checkFileExtension(it, extension) {
  // Check file extension
  if (!it.toString().toLowerCase().endsWith(extension.toLowerCase())) exit 1, "File: ${it} has the wrong extension: ${extension} see --help for more information"
}

def checkParameterExistence(it, list) {
  // Check parameter existence
  if (!list.contains(it)) {
    println("Unknown parameter: ${it}")
    return false
  }
  return true
}

def checkParamReturnFile(item) {
  params."${item}" = params.genomes[params.genome]."${item}"
  return file(params."${item}")
}

def checkReferenceMap(referenceMap) {
  // Loop through all the references files to check their existence
  referenceMap.every {
    referenceFile, fileToCheck ->
    checkRefExistence(referenceFile, fileToCheck)
  }
}

def checkRefExistence(referenceFile, fileToCheck) {
  if (fileToCheck instanceof List) return fileToCheck.every{ checkRefExistence(referenceFile, it) }
  def f = file(fileToCheck)
  // this is an expanded wildcard: we can assume all files exist
  if (f instanceof List && f.size() > 0) return true
  else if (!f.exists()) {
    log.info  "Missing references: ${referenceFile} ${fileToCheck}"
    return false
  }
  return true
}

def checkUppmaxProject() {
  // check if UPPMAX project number is specified
  return !(workflow.profile == 'slurm' && !params.project)
}

def checkExactlyOne(list) {
  final n = 0
  list.each{n += it ? 1 : 0}
  return n == 1
}

def defineDirectoryMap() {
  return [
    'nonRealigned'     : "${params.outDir}/Preprocessing/NonRealigned",
    'nonRecalibrated'  : "${params.outDir}/Preprocessing/NonRecalibrated",
    'recalibrated'     : "${params.outDir}/Preprocessing/Recalibrated",
    'bamQC'            : "${params.outDir}/Reports/bamQC",
    'bcftoolsStats'    : "${params.outDir}/Reports/BCFToolsStats",
    'fastQC'           : "${params.outDir}/Reports/FastQC",
    'markDuplicatesQC' : "${params.outDir}/Reports/MarkDuplicates",
    'samtoolsStats'    : "${params.outDir}/Reports/SamToolsStats"
  ]
}

def defineReferenceMap() {
  if (!(params.genome in params.genomes)) exit 1, "Genome ${params.genome} not found in configuration"
  return [
    'dbsnp'            : checkParamReturnFile("dbsnp"),
    'dbsnpIndex'       : checkParamReturnFile("dbsnpIndex"),
    // genome reference dictionary
    'genomeDict'       : checkParamReturnFile("genomeDict"),
    // FASTA genome reference
    'genomeFile'       : checkParamReturnFile("genomeFile"),
    // genome .fai file
    'genomeIndex'      : checkParamReturnFile("genomeIndex"),
    // BWA index files
    'bwaIndex'         : checkParamReturnFile("bwaIndex"),
    // intervals file for spread-and-gather processes
    'intervals'        : checkParamReturnFile("intervals"),
    // VCFs with known indels (such as 1000 Genomes, Mill’s gold standard)
    'knownIndels'      : checkParamReturnFile("knownIndels"),
    'knownIndelsIndex' : checkParamReturnFile("knownIndelsIndex"),
  ]
}

def defineStepList() {
  return [
    'mapping',
    'realign',
    'recalibrate'
  ]
}

def extractBams(tsvFile) {
  // Channeling the TSV file containing BAM.
  // Format is: "subject gender status sample bam bai"
  Channel
    .from(tsvFile.readLines())
    .map{line ->
      def list      = returnTSV(line.split(),6)
      def idPatient = list[0]
      def gender    = list[1]
      def status    = returnStatus(list[2].toInteger())
      def idSample  = list[3]
      def bamFile   = returnFile(list[4])
      def baiFile   = returnFile(list[5])

      checkFileExtension(bamFile,".bam")
      checkFileExtension(baiFile,".bai")

      [ idPatient, gender, status, idSample, bamFile, baiFile ]
    }
}

def extractFastq(tsvFile) {
  // Channeling the TSV file containing FASTQ.
  // Format is: "subject gender status sample lane fastq1 fastq2"
  Channel
    .from(tsvFile.readLines())
    .map{line ->
      def list       = returnTSV(line.split(),7)
      def idPatient  = list[0]
      def gender     = list[1]
      def status     = returnStatus(list[2].toInteger())
      def idSample   = list[3]
      def idRun      = list[4]

      // Normally path to files starts from workflow.launchDir
      // But when executing workflow from Github
      // Path to hosted FASTQ files starts from workflow.projectDir
      def fastqFile1 = workflow.commitId && params.test ? returnFile("${workflow.projectDir}/${list[5]}") : returnFile(list[5])
      def fastqFile2 = workflow.commitId && params.test ? returnFile("${workflow.projectDir}/${list[6]}") : returnFile(list[6])

      checkFileExtension(fastqFile1,".fastq.gz")
      checkFileExtension(fastqFile2,".fastq.gz")

      [idPatient, gender, status, idSample, idRun, fastqFile1, fastqFile2]
    }
}

def extractFastqFromDir(pattern) {
  // create a channel of FASTQs from a directory pattern such as
  // "my_samples/*/". All samples are considered 'normal'.
  // All FASTQ files in subdirectories are collected and emitted;
  // they must have _R1_ and _R2_ in their names.

  def fastq = Channel.create()

  // a temporary channel does all the work
  Channel
    .fromPath(pattern, type: 'dir')
    .ifEmpty { error "No directories found matching pattern '${pattern}'" }
    .subscribe onNext: { sampleDir ->
      // the last name of the sampleDir is assumed to be a unique sample id
      sampleId = sampleDir.getFileName().toString()

      for (path1 in file("${sampleDir}/**_R1_*.fastq.gz")) {
        assert path1.getName().contains('_R1_')
        path2 = file(path1.toString().replace('_R1_', '_R2_'))
        if (!path2.exists()) error "Path '${path2}' not found"
        (flowcell, lane) = flowcellLaneFromFastq(path1)
        patient = sampleId
        gender = 'ZZ'  // unused
        status = 0  // normal (not tumor)
        rgId = "${flowcell}.${sampleId}.${lane}"
        result = [patient, gender, status, sampleId, rgId, path1, path2]
        fastq.bind(result)
      }
  }, onComplete: { fastq.close() }

  fastq
}

def extractRecal(tsvFile) {
  // Channeling the TSV file containing Recalibration Tables.
  // Format is: "subject gender status sample bam bai recalTables"
  Channel
    .from(tsvFile.readLines())
    .map{line ->
      def list       = returnTSV(line.split(),7)
      def idPatient  = list[0]
      def gender     = list[1]
      def status     = returnStatus(list[2].toInteger())
      def idSample   = list[3]
      def bamFile    = returnFile(list[4])
      def baiFile    = returnFile(list[5])
      def recalTable = returnFile(list[6])

      checkFileExtension(bamFile,".bam")
      checkFileExtension(baiFile,".bai")
      checkFileExtension(recalTable,".recal.table")

      [ idPatient, gender, status, idSample, bamFile, baiFile, recalTable ]
    }
}

def extractGenders(channel) {
  def genders = [:]  // an empty map
  channel = channel.map{ it ->
    def idPatient = it[0]
    def gender = it[1]
    genders[idPatient] = gender
    [idPatient] + it[2..-1]
  }
  [genders, channel]
}

def flowcellLaneFromFastq(path) {
  // parse first line of a FASTQ file (optionally gzip-compressed)
  // and return the flowcell id and lane number.
  // expected format:
  // xx:yy:FLOWCELLID:LANE:... (seven fields)
  // or
  // FLOWCELLID:LANE:xx:... (five fields)
  InputStream fileStream = new FileInputStream(path.toFile())
  InputStream gzipStream = new java.util.zip.GZIPInputStream(fileStream)
  Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
  BufferedReader buffered = new BufferedReader(decoder)
  def line = buffered.readLine()
  assert line.startsWith('@')
  line = line.substring(1)
  def fields = line.split(' ')[0].split(':')
  String fcid
  int lane
  if (fields.size() == 7) {
    // CASAVA 1.8+ format
    fcid = fields[2]
    lane = fields[3].toInteger()
  }
  else if (fields.size() == 5) {
    fcid = fields[0]
    lane = fields[1].toInteger()
  }
  [fcid, lane]
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
}

def minimalInformationMessage() {
  // Minimal information message
  log.info "Command Line: " + workflow.commandLine
  log.info "Profile     : " + workflow.profile
  log.info "Project Dir : " + workflow.projectDir
  log.info "Launch Dir  : " + workflow.launchDir
  log.info "Work Dir    : " + workflow.workDir
  log.info "Out Dir     : " + params.outDir
  log.info "TSV file    : ${tsvFile}"
  log.info "Genome      : " + params.genome
  log.info "Genome_base : " + params.genome_base
  log.info "Step        : " + step
  log.info "Containers  :"
  if (params.repository) log.info "  Repository   : ${params.repository}"
  else log.info "  ContainerPath: " + params.containerPath
  log.info "  Tag          : " + params.tag
  log.info "Reference files used:"
  log.info "  dbsnp       :\n\t" + referenceMap.dbsnp
  log.info "\t" + referenceMap.dbsnpIndex
  log.info "  genome      :\n\t" + referenceMap.genomeFile
  log.info "\t" + referenceMap.genomeDict
  log.info "\t" + referenceMap.genomeIndex
  log.info "  bwa indexes :\n\t" + referenceMap.bwaIndex.join(',\n\t')
  log.info "  intervals   :\n\t" + referenceMap.intervals
  log.info "  knownIndels :\n\t" + referenceMap.knownIndels.join(',\n\t')
  log.info "\t" + referenceMap.knownIndelsIndex.join(',\n\t')
}

def nextflowMessage() {
  // Nextflow message (version + build)
  log.info "N E X T F L O W  ~  version ${workflow.nextflow.version} ${workflow.nextflow.build}"
}

def returnFile(it) {
  // return file if it exists
  if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
  return file(it)
}

def returnStatus(it) {
  // Return status if it's correct
  // Status should be only 0 or 1
  // 0 being normal
  // 1 being tumor (or relapse or anything that is not normal...)
  if (!(it in [0, 1])) exit 1, "Status is not recognized in TSV file: ${it}, see --help for more information"
  return it
}

def returnTSV(it, number) {
  // return TSV if it has the correct number of items in row
  if (it.size() != number) exit 1, "Malformed row in TSV file: ${it}, see --help for more information"
  return it
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
