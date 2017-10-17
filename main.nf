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
 - CreateIntervalBeds - Create and sort intervals into bed files
 - RunHaplotypecaller - Run HaplotypeCaller for Germline Variant Calling (Parallelized processes)
 - RunGenotypeGVCFs - Run HaplotypeCaller for Germline Variant Calling (Parallelized processes)
 - RunMutect1 - Run MuTect1 for Variant Calling (Parallelized processes)
 - RunMutect2 - Run MuTect2 for Variant Calling (Parallelized processes)
 - RunFreeBayes - Run FreeBayes for Variant Calling (Parallelized processes)
 - ConcatVCF - Merge results from HaplotypeCaller, MuTect1 and MuTect2
 - RunStrelka - Run Strelka for Variant Calling
 - RunSingleStrelka - Run Strelka for Germline Variant Calling
 - RunManta - Run Manta for Structural Variant Calling
 - RunSingleManta - Run Manta for Single Structural Variant Calling
 - RunAlleleCount - Run AlleleCount to prepare for ASCAT
 - RunConvertAlleleCounts - Run convertAlleleCounts to prepare for ASCAT
 - RunAscat - Run ASCAT for CNV
 - RunBcftoolsStats - Run BCFTools stats on vcf before annotation
 - RunSnpeff - Run snpEff for annotation of vcf files
 - RunVEP - Run VEP for annotation of vcf files
 - GenerateMultiQCconfig - Generate a config file for MultiQC
 - RunMultiQC - Run MultiQC for report and QC
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
    log.error "============================================================\n" +
              "  Nextflow version $nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please update Nextflow.\n" +
              "============================================================"
}

// Default params:
// Such params are overridden by command line or configuration definitions

// No tools to annotate
params.annotateTools = ''
// No vcf to annotare
params.annotateVCF = ''
// For MultiQC reports
params.callName = ''
// For MultiQC reports
params.contactMail = ''
// GVCF are generated
params.noGVCF = false
// Reports are generated
params.noReports = false
// No sample is defined
params.sample = ''
// No sampleDir is defined
params.sampleDir = ''
// Step is mapping
params.step = 'mapping'
// Not testing
params.test = ''
// No tools to be used
params.tools = ''

if (params.help) exit 0, helpMessage()
if (params.version) exit 0, versionMessage()
if (!isAllowedParams(params)) exit 1, "params is unknown, see --help for more information"
if (!checkUppmaxProject()) exit 1, "No UPPMAX project ID found! Use --project <UPPMAX Project ID>"

step = params.step.toLowerCase()
if (step == 'preprocessing') step = 'mapping'
if (step == 'skippreprocessing') step = 'variantcalling'
tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase()} : []
annotateTools = params.annotateTools ? params.annotateTools.split(',').collect{it.trim().toLowerCase()} : []
annotateVCF = params.annotateVCF ? params.annotateVCF.split(',').collect{it.trim()} : []

directoryMap = defineDirectoryMap()
referenceMap = defineReferenceMap()
stepList = defineStepList()
toolList = defineToolList()
nucleotidesPerSecond = 1000.0 // used to estimate variant calling runtime
gvcf = !params.noGVCF
reports = !params.noReports
verbose = params.verbose

if (!checkParameterExistence(step, stepList)) exit 1, 'Unknown step, see --help for more information'
if (step.contains(',')) exit 1, 'You can choose only one step, see --help for more information'
if (step == 'mapping' && !checkExactlyOne([params.test, params.sample, params.sampleDir]))
  exit 1, 'Please define which samples to work on by providing exactly one of the --test, --sample or --sampleDir options'
if (!checkReferenceMap(referenceMap)) exit 1, 'Missing Reference file(s), see --help for more information'
if (!checkParameterList(tools,toolList)) exit 1, 'Unknown tool(s), see --help for more information'

if (params.test && params.genome in ['GRCh37', 'GRCh38']) {
  referenceMap.intervals = file("$workflow.projectDir/repeats/tiny_${params.genome}.list")
}

// TODO
// MuTect and Mutect2 could be run without a recalibrated BAM (they support
// the --BQSR option), but this is not implemented, yet.
// TODO
// FreeBayes does not need recalibrated BAMs, but we need to test whether
// the channels are set up correctly when we disable it
if (step == "recalibrate" && tools != ['haplotypecaller']) explicitBqsrNeeded = true
else explicitBqsrNeeded = tools.intersect(['manta', 'mutect1', 'mutect2', 'vardict', 'freebayes', 'strelka']).asBoolean()

tsvPath = ''
if (params.sample) tsvPath = params.sample

// No need for tsv file for step annotate
if (!params.sample && !params.sampleDir) {
  tsvPaths = [
  'mapping': "$workflow.projectDir/data/tsv/tiny.tsv",
  'realign': "$workflow.launchDir/$directoryMap.nonRealigned/nonRealigned.tsv",
  'recalibrate': "$workflow.launchDir/$directoryMap.nonRecalibrated/nonRecalibrated.tsv",
  'variantcalling': "$workflow.launchDir/$directoryMap.recalibrated/recalibrated.tsv"
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
    case 'variantcalling': bamFiles = extractBams(tsvFile); break
    default: exit 1, "Unknown step $step"
  }
} else if (params.sampleDir) {
  if (step != 'mapping') exit 1, '--sampleDir does not support steps other than "mapping"'
  fastqFiles = extractFastqFromDir(params.sampleDir)
  tsvFile = params.sampleDir  // used in the reports
} else if (step != 'annotate') exit 1, 'No sample were defined, see --help'

if (step == 'mapping') {
  (patientGenders, fastqFiles) = extractGenders(fastqFiles)
} else {
  (patientGenders, bamFiles) = extractGenders(bamFiles)
}

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

startMessage()

(fastqFiles, fastqFilesforFastQC) = fastqFiles.into(2)

if (verbose) fastqFiles = fastqFiles.view {
  "FASTQs to preprocess:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\tRun   : ${it[3]}\n\
  Files : [${it[4].fileName}, ${it[5].fileName}]"
}

if (verbose) bamFiles = bamFiles.view {
  "BAMs to process:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  Files : [${it[3].fileName}, ${it[4].fileName}]"
}

process RunFastQC {
  tag {idPatient + "-" + idRun}

  publishDir "$directoryMap.fastQC/$idRun", mode: 'copy'

  input:
    set idPatient, status, idSample, idRun, file(fastqFile1), file(fastqFile2) from fastqFilesforFastQC

  output:
    file "*_fastqc.{zip,html}" into fastQCreport

  when: step == 'mapping' && reports

  script:
  """
  fastqc -q $fastqFile1 $fastqFile2
  """
}

if (verbose) fastQCreport = fastQCreport.view {
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

  when: step == 'mapping'

  script:
  readGroup = "@RG\\tID:$idRun\\tPU:$idRun\\tSM:$idSample\\tLB:$idSample\\tPL:illumina"
  // adjust mismatch penalty for tumor samples
  extra = status == 1 ? "-B 3 " : ""
  """
  bwa mem -R \"$readGroup\" ${extra}-t $task.cpus -M \
  $genomeFile $fastqFile1 $fastqFile2 | \
  samtools sort --threads $task.cpus -m 4G - > ${idRun}.bam
  """
}

if (verbose) mappedBam = mappedBam.view {
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

  when: step == 'mapping'

  script:
  """
  samtools merge --threads $task.cpus ${idSample}.bam $bam
  """
}

if (verbose) singleBam = singleBam.view {
  "Single BAM:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  File  : [${it[3].fileName}]"
}

if (verbose) mergedBam = mergedBam.view {
  "Merged BAM:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  File  : [${it[3].fileName}]"
}

mergedBam = mergedBam.mix(singleBam)

if (verbose) mergedBam = mergedBam.view {
  "BAM for MarkDuplicates:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  File  : [${it[3].fileName}]"
}

process MarkDuplicates {
  tag {idPatient + "-" + idSample}

  publishDir '.', saveAs: { it == "${bam}.metrics" ? "$directoryMap.markDuplicatesQC/$it" : "$directoryMap.nonRealigned/$it" }, mode: 'copy'

  input:
    set idPatient, status, idSample, file(bam) from mergedBam

  output:
    set idPatient, file("${idSample}_${status}.md.bam"), file("${idSample}_${status}.md.bai") into duplicates
    set idPatient, status, idSample, val("${idSample}_${status}.md.bam"), val("${idSample}_${status}.md.bai") into markDuplicatesTSV
    file ("${bam}.metrics") into markDuplicatesReport

  when: step == 'mapping'

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
  "$idPatient\t$gender\t$status\t$idSample\t$directoryMap.nonRealigned/$bam\t$directoryMap.nonRealigned/$bai\n"
}.collectFile(
  name: 'nonRealigned.tsv', sort: true, storeDir: directoryMap.nonRealigned
)

// Create intervals for realignement using both tumor+normal as input
// Group the marked duplicates BAMs for intervals and realign by idPatient
// Grouping also by gender, to make a nicer channel
duplicatesGrouped = Channel.empty()
if (step == 'mapping') {
  duplicatesGrouped = duplicates.groupTuple()
} else if (step == 'realign') {
  duplicatesGrouped = bamFiles.map{
    idPatient, status, idSample, bam, bai ->
    [idPatient, bam, bai]
  }.groupTuple()
}

// The duplicatesGrouped channel is duplicated
// one copy goes to the RealignerTargetCreator process
// and the other to the IndelRealigner process
(duplicatesInterval, duplicatesRealign) = duplicatesGrouped.into(2)

if (verbose) duplicatesInterval = duplicatesInterval.view {
  "BAMs for RealignerTargetCreator:\n\
  ID    : ${it[0]}\n\
  Files : ${it[1].fileName}\n\
  Files : ${it[2].fileName}"
}

if (verbose) duplicatesRealign = duplicatesRealign.view {
  "BAMs to phase:\n\
  ID    : ${it[0]}\n\
  Files : ${it[1].fileName}\n\
  Files : ${it[2].fileName}"
}

if (verbose) markDuplicatesReport = markDuplicatesReport.view {
  "MarkDuplicates report:\n\
  File  : [$it.fileName]"
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

  when: step == 'mapping' || step == 'realign'

  script:
  bams = bam.collect{"-I $it"}.join(' ')
  known = knownIndels.collect{"-known $it"}.join(' ')
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar \$GATK_HOME/GenomeAnalysisTK.jar \
  -T RealignerTargetCreator \
  $bams \
  -R $genomeFile \
  $known \
  -nt $task.cpus \
  -L $intervals \
  -o ${idPatient}.intervals
  """
}

if (verbose) intervals = intervals.view {
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

if (verbose) bamsAndIntervals = bamsAndIntervals.view {
  "BAMs and Intervals phased for IndelRealigner:\n\
  ID    : ${it[0]}\n\
  Files : ${it[1].fileName}\n\
  Files : ${it[2].fileName}\n\
  File  : [${it[3].fileName}]"
}

// use nWayOut to split into T/N pair again
process IndelRealigner {
  tag {idPatient}

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

  when: step == 'mapping' || step == 'realign'

  script:
  bams = bam.collect{"-I $it"}.join(' ')
  known = knownIndels.collect{"-known $it"}.join(' ')
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar \$GATK_HOME/GenomeAnalysisTK.jar \
  -T IndelRealigner \
  $bams \
  -R $genomeFile \
  -targetIntervals $intervals \
  $known \
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

if (verbose) realignedBam = realignedBam.view {
  "Realigned BAM to CreateRecalibrationTable:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  Files : [${it[3].fileName}, ${it[4].fileName}]"
}

process CreateRecalibrationTable {
  tag {idPatient + "-" + idSample}

  publishDir directoryMap.nonRecalibrated, mode: 'copy'

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

  when: step == 'mapping' || step == 'realign'

  script:
  known = knownIndels.collect{ "-knownSites $it" }.join(' ')
  """
  java -Xmx${task.memory.toGiga()}g \
  -Djava.io.tmpdir="/tmp" \
  -jar \$GATK_HOME/GenomeAnalysisTK.jar \
  -T BaseRecalibrator \
  -R $genomeFile \
  -I $bam \
  -L $intervals \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  -knownSites $dbsnp \
  $known \
  -nct $task.cpus \
  -l INFO \
  -o ${idSample}.recal.table
  """
}
// Create a TSV file to restart from this step
recalibrationTableTSV.map { idPatient, status, idSample, bam, bai, recalTable ->
  gender = patientGenders[idPatient]
  "$idPatient\t$gender\t$status\t$idSample\t$directoryMap.nonRecalibrated/$bam\t$directoryMap.nonRecalibrated/$bai\t\t$directoryMap.nonRecalibrated/$recalTable\n"
}.collectFile(
  name: 'nonRecalibrated.tsv', sort: true, storeDir: directoryMap.nonRecalibrated
)

if (step == 'recalibrate') recalibrationTable = bamFiles

if (verbose) recalibrationTable = recalibrationTable.view {
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

  publishDir directoryMap.recalibrated, mode: 'copy'

  input:
    set idPatient, status, idSample, file(bam), file(bai), recalibrationReport from recalibrationTable
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
  when: step != 'variantcalling' && explicitBqsrNeeded

  script:
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar \$GATK_HOME/GenomeAnalysisTK.jar \
  -T PrintReads \
  -R $genomeFile \
  -I $bam \
  -L $intervals \
  --BQSR $recalibrationReport \
  -o ${idSample}.recal.bam
  """
}
// Creating a TSV file to restart from this step
recalibratedBamTSV.map { idPatient, status, idSample, bam, bai ->
  gender = patientGenders[idPatient]
  "$idPatient\t$gender\t$status\t$idSample\t$directoryMap.recalibrated/$bam\t$directoryMap.recalibrated/$bai\n"
}.collectFile(
  name: 'recalibrated.tsv', sort: true, storeDir: directoryMap.recalibrated
)

if (step == 'variantcalling') {
  // assume input is recalibrated, ignore explicitBqsrNeeded
  (recalibratedBam, recalTables) = bamFiles.into(2)

  recalTables = recalTables.map{ it + [null] } // null recalibration table means: do not use --BQSR

  (recalTables, recalibrationTableForHC) = recalTables.into(2)
  recalTables = recalTables.map { [it[0]] + it[2..-1] } // remove status
} else if (!explicitBqsrNeeded) {
  (bamForBamQC, bamForSamToolsStats, recalibratedBam) = recalibrationTableForHC.map { it[0..-2] }.into(3)
}

if (verbose) recalibratedBam = recalibratedBam.view {
  "Recalibrated BAM for variant Calling:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  Files : [${it[3].fileName}, ${it[4].fileName}]"
}

process RunSamtoolsStats {
  tag {idPatient + "-" + idSample}

  publishDir directoryMap.samtoolsStats, mode: 'copy'

  input:
    set idPatient, status, idSample, file(bam), file(bai) from bamForSamToolsStats

  output:
    file ("${bam}.samtools.stats.out") into samtoolsStatsReport

    when: reports

    script:
    """
    samtools stats $bam > ${bam}.samtools.stats.out
    """
}

if (verbose) samtoolsStatsReport = samtoolsStatsReport.view {
  "SAMTools stats report:\n\
  File  : [${it.fileName}]"
}

process RunBamQC {
  tag {idPatient + "-" + idSample}

  publishDir directoryMap.bamQC, mode: 'copy'

  input:
    set idPatient, status, idSample, file(bam), file(bai) from bamForBamQC

  output:
    file("$idSample") into bamQCreport

    when: reports

    script:
    """
    qualimap --java-mem-size=${task.memory.toGiga()}G bamqc -bam $bam -outdir $idSample -outformat HTML
    """
}

if (verbose) bamQCreport = bamQCreport.view {
  "BamQC report:\n\
  Dir   : [${it.fileName}]"
}

// Here we have a recalibrated bam set, but we need to separate the bam files based on patient status.
// The sample tsv config file which is formatted like: "subject status sample lane fastq1 fastq2"
// cf fastqFiles channel, I decided just to add _status to the sample name to have less changes to do.
// And so I'm sorting the channel if the sample match _0, then it's a normal sample, otherwise tumor.
// Then combine normal and tumor to get each possibilities
// ie. normal vs tumor1, normal vs tumor2, normal vs tumor3
// then copy this channel into channels for each variant calling
// I guess it will still work even if we have multiple normal samples

// separate recalibrateBams by status
bamsNormal = Channel.create()
bamsTumor = Channel.create()

recalibratedBam
  .choice(bamsTumor, bamsNormal) {it[1] == 0 ? 1 : 0}

// Ascat, Strelka Germline & Manta Germline SV
bamsForAscat = Channel.create()
bamsForSingleManta = Channel.create()
bamsForSingleStrelka = Channel.create()

(bamsTumorTemp, bamsTumor) = bamsTumor.into(2)
(bamsNormalTemp, bamsNormal) = bamsNormal.into(2)
(bamsForAscat, bamsForSingleManta, bamsForSingleStrelka) = bamsNormalTemp.mix(bamsTumorTemp).into(3)

// Removing status because not relevant anymore
bamsNormal = bamsNormal.map { idPatient, status, idSample, bam, bai -> [idPatient, idSample, bam, bai] }

bamsTumor = bamsTumor.map { idPatient, status, idSample, bam, bai -> [idPatient, idSample, bam, bai] }

// We know that MuTect2 (and other somatic callers) are notoriously slow.
// To speed them up we are chopping the reference into smaller pieces.
// (see repeats/centromeres.list).
// Do variant calling by this intervals, and re-merge the VCFs.
// Since we are on a cluster, this can parallelize the variant call processes.
// And push down the variant call wall clock time significanlty.

process CreateIntervalBeds {
  tag {intervals.fileName}

  input:
    file(intervals) from Channel.value(referenceMap.intervals)

  output:
    file '*.bed' into bedIntervals mode flatten

  script:
  // If the interval file is BED format, the fifth column is interpreted to
  // contain runtime estimates, which is then used to combine short-running jobs
  if (intervals.getName().endsWith('.bed'))
    """
    awk -vFS="\t" '{
      t = \$5  # runtime estimate
      if (t == "") {
        # no runtime estimate in this row, assume default value
        t = (\$3 - \$2) / ${nucleotidesPerSecond}
      }
      if (name == "" || (chunk > 600 && (chunk + t) > longest * 1.05)) {
        # start a new chunk
        name = sprintf("%s_%d-%d.bed", \$1, \$2+1, \$3)
        chunk = 0
        longest = 0
      }
      if (t > longest)
        longest = t
      chunk += t
      print \$0 > name
    }' $intervals
    """
  else
    """
    awk -vFS="[:-]" '{
      name = sprintf("%s_%d-%d", \$1, \$2, \$3);
      printf("%s\\t%d\\t%d\\n", \$1, \$2-1, \$3) > name ".bed"
    }' $intervals
    """
}

bedIntervals = bedIntervals
  .map { intervalFile ->
    final duration = 0.0
    for (line in intervalFile.readLines()) {
      final fields = line.split('\t')
      if (fields.size() >= 5) {
        duration += fields[4].toFloat()
      } else {
        start = fields[1].toInteger()
        end = fields[2].toInteger()
        duration += (end - start) / nucleotidesPerSecond
      }
    }
    [duration, intervalFile]
  }.toSortedList({ a, b -> b[0] <=> a[0] })
  .flatten().collate(2)
  .map{duration, intervalFile -> intervalFile}

if (verbose) bedIntervals = bedIntervals.view {
  "  Interv: ${it.baseName}"
}

(bamsNormalTemp, bamsNormal, bedIntervals) = generateIntervalsForVC(bamsNormal, bedIntervals)
(bamsTumorTemp, bamsTumor, bedIntervals) = generateIntervalsForVC(bamsTumor, bedIntervals)

// HaplotypeCaller
bamsForHC = bamsNormalTemp.mix(bamsTumorTemp)
bedIntervals = bedIntervals.tap { intervalsTemp }
recalTables = recalTables
  .spread(intervalsTemp)
  .map { patient, sample, bam, bai, recalTable, intervalBed ->
    [patient, sample, bam, bai, intervalBed, recalTable] }

// re-associate the BAMs and samples with the recalibration table
bamsForHC = bamsForHC
  .phase(recalTables) { it[0..4] }
  .map { it1, it2 -> it1 + [it2[6]] }

bamsAll = bamsNormal.combine(bamsTumor)

// Since idPatientNormal and idPatientTumor are the same
// It's removed from bamsAll Channel (same for genderNormal)
// /!\ It is assumed that every sample are from the same patient
bamsAll = bamsAll.map {
  idPatientNormal, idSampleNormal, bamNormal, baiNormal, idPatientTumor, idSampleTumor, bamTumor, baiTumor ->
  [idPatientNormal, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor]
}

// Manta and Strelka
(bamsForManta, bamsForStrelka, bamsAll) = bamsAll.into(3)

bamsTumorNormalIntervals = bamsAll.spread(bedIntervals)

// MuTect1, MuTect2, FreeBayes
(bamsFMT1, bamsFMT2, bamsFFB) = bamsTumorNormalIntervals.into(3)


process RunHaplotypecaller {
  tag {idSample + "-" + intervalBed.baseName}

  input:
    set idPatient, idSample, file(bam), file(bai), file(intervalBed), recalTable from bamsForHC //Are these values `ped to bamNormal already?
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(dbsnp), file(dbsnpIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.dbsnp,
      referenceMap.dbsnpIndex
    ])

  output:
    set val("gvcf-hc"), idPatient, idSample, idSample, file("${intervalBed.baseName}_${idSample}.g.vcf") into hcGenomicVCF
    set idPatient, idSample, file(intervalBed), file("${intervalBed.baseName}_${idSample}.g.vcf") into vcfsToGenotype

  when: 'haplotypecaller' in tools

  script:
  BQSR = (recalTable != null) ? "--BQSR $recalTable" : ''
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar \$GATK_HOME/GenomeAnalysisTK.jar \
  -T HaplotypeCaller \
  --emitRefConfidence GVCF \
  -pairHMM LOGLESS_CACHING \
  -R $genomeFile \
  --dbsnp $dbsnp \
  $BQSR \
  -I $bam \
  -L $intervalBed \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  -o ${intervalBed.baseName}_${idSample}.g.vcf
  """
}
hcGenomicVCF = hcGenomicVCF.groupTuple(by:[0,1,2,3])

if (!gvcf) {hcGenomicVCF.close()}

process RunGenotypeGVCFs {
  tag {idSample + "-" + intervalBed.baseName}

  input:
    set idPatient, idSample, file(intervalBed), file(gvcf) from vcfsToGenotype
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(dbsnp), file(dbsnpIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.dbsnp,
      referenceMap.dbsnpIndex
    ])

  output:
    set val("haplotypecaller"), idPatient, idSample, idSample, file("${intervalBed.baseName}_${idSample}.vcf") into hcGenotypedVCF

  when: 'haplotypecaller' in tools

  script:
  // Using -L is important for speed
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar \$GATK_HOME/GenomeAnalysisTK.jar \
  -T GenotypeGVCFs \
  -R $genomeFile \
  -L $intervalBed \
  --dbsnp $dbsnp \
  --variant $gvcf \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  -o ${intervalBed.baseName}_${idSample}.vcf
  """
}
hcGenotypedVCF = hcGenotypedVCF.groupTuple(by:[0,1,2,3])

process RunMutect1 {
  tag {idSampleTumor + "_vs_" + idSampleNormal + "-" + intervalBed.baseName}

  input:
    set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), file(intervalBed) from bamsFMT1
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(dbsnp), file(dbsnpIndex), file(cosmic), file(cosmicIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.dbsnp,
      referenceMap.dbsnpIndex,
      referenceMap.cosmic,
      referenceMap.cosmicIndex
    ])

  output:
    set val("mutect1"), idPatient, idSampleNormal, idSampleTumor, file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf") into mutect1Output

  when: 'mutect1' in tools

  script:
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar \$MUTECT_HOME/muTect.jar \
  -T MuTect \
  -R $genomeFile \
  --cosmic $cosmic \
  --dbsnp $dbsnp \
  -I:normal $bamNormal \
  -I:tumor $bamTumor \
  -L $intervalBed \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --out ${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.call_stats.out \
  --vcf ${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf
  """
}

mutect1Output = mutect1Output.groupTuple(by:[0,1,2,3])

process RunMutect2 {
  tag {idSampleTumor + "_vs_" + idSampleNormal + "-" + intervalBed.baseName}

  input:
    set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), file(intervalBed) from bamsFMT2
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(dbsnp), file(dbsnpIndex), file(cosmic), file(cosmicIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.dbsnp,
      referenceMap.dbsnpIndex,
      referenceMap.cosmic,
      referenceMap.cosmicIndex
    ])

  output:
    set val("mutect2"), idPatient, idSampleNormal, idSampleTumor, file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf") into mutect2Output

  when: 'mutect2' in tools

  script:
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar \$GATK_HOME/GenomeAnalysisTK.jar \
  -T MuTect2 \
  -R $genomeFile \
  --cosmic $cosmic \
  --dbsnp $dbsnp \
  -I:normal $bamNormal \
  -I:tumor $bamTumor \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  -L $intervalBed \
  -o ${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf
  """
}

mutect2Output = mutect2Output.groupTuple(by:[0,1,2,3])

process RunFreeBayes {
  tag {idSampleTumor + "_vs_" + idSampleNormal + "-" + intervalBed.baseName}

  input:
    set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), file(intervalBed) from bamsFFB
    file(genomeFile) from Channel.value(referenceMap.genomeFile)

  output:
    set val("freebayes"), idPatient, idSampleNormal, idSampleTumor, file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf") into freebayesOutput

  when: 'freebayes' in tools

  script:
  """
  freebayes \
    -f $genomeFile \
    --pooled-continuous \
    --pooled-discrete \
    --genotype-qualities \
    --report-genotype-likelihood-max \
    --allele-balance-priors-off \
    --min-alternate-fraction 0.03 \
    --min-repeat-entropy 1 \
    --min-alternate-count 2 \
    -t $intervalBed \
    $bamTumor \
    $bamNormal > ${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf
  """
}

freebayesOutput = freebayesOutput.groupTuple(by:[0,1,2,3])

// we are merging the VCFs that are called separatelly for different intervals
// so we can have a single sorted VCF containing all the calls for a given caller

vcfsToMerge = hcGenomicVCF.mix(hcGenotypedVCF, mutect1Output, mutect2Output, freebayesOutput)
if (verbose) vcfsToMerge = vcfsToMerge.view {
  "VCFs To be merged:\n\
  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: [${it[3]}, ${it[2]}]\n\
  Files : ${it[4].fileName}"
}

process ConcatVCF {
  tag {variantCaller in ['gvcf-hc', 'haplotypecaller'] ? variantCaller + "-" + idSampleNormal : variantCaller + "_" + idSampleTumor + "_vs_" + idSampleNormal}

  publishDir "${directoryMap."$variantCaller"}", mode: 'copy'

  input:
    set variantCaller, idPatient, idSampleNormal, idSampleTumor, file(vcFiles) from vcfsToMerge
    file(genomeIndex) from Channel.value(referenceMap.genomeIndex)

  output:
    set variantCaller, idPatient, idSampleNormal, idSampleTumor, file("*.vcf.gz") into vcfConcatenated
    file("*.vcf.gz.tbi") into vcfConcatenatedTbi

  when: 'haplotypecaller' in tools || 'mutect1' in tools || 'mutect2' in tools || 'freebayes' in tools

  script:
  if (variantCaller == 'haplotypecaller') {
    outputFile = "${variantCaller}_${idSampleNormal}.vcf"
  } else if (variantCaller == 'gvcf-hc') {
    outputFile = "haplotypecaller_${idSampleNormal}.g.vcf"
  } else {
    outputFile = "${variantCaller}_${idSampleTumor}_vs_${idSampleNormal}.vcf"
  }
  vcfFiles = vcFiles.collect{" $it"}.join(' ')

  """
  # first make a header from one of the VCF intervals
  # get rid of interval information only from the GATK command-line, but leave the rest
  FIRSTVCF=\$(ls *.vcf | head -n 1)
  sed -n '/^[^#]/q;p' \$FIRSTVCF | \
  awk '!/GATKCommandLine/{print}/GATKCommandLine/{for(i=1;i<=NF;i++){if(\$i!~/intervals=/ && \$i !~ /out=/){printf("%s ",\$i)}}printf("\\n")}' \
  > header

  # Get list of contigs from the FASTA index (.fai). We cannot use the ##contig
  # header in the VCF as it is optional (FreeBayes does not save it, for example)
  CONTIGS=(\$(cut -f1 ${genomeIndex}))

  # concatenate VCFs in the correct order
  (
    cat header

    for chr in "\${CONTIGS[@]}"; do
      # Skip if globbing would not match any file to avoid errors such as
      # "ls: cannot access chr3_*.vcf: No such file or directory" when chr3
      # was not processed.
      pattern="\${chr}_*.vcf"
      if ! compgen -G "\${pattern}" > /dev/null; then continue; fi

      # ls -v sorts by numeric value ("version"), which means that chr1_100_
      # is sorted *after* chr1_99_.
      for vcf in \$(ls -v \${pattern}); do
        # Determine length of header.
        # The 'q' command makes sed exit when it sees the first non-header
        # line, which avoids reading in the entire file.
        L=\$(sed -n '/^[^#]/q;p' \${vcf} | wc -l)

        # Then print all non-header lines. Since tail is very fast (nearly as
        # fast as cat), this is way more efficient than using a single sed,
        # awk or grep command.
        tail -n +\$((L+1)) \${vcf}
      done
    done
  ) | bgzip > ${outputFile}.gz
  tabix ${outputFile}.gz
  """
}

if (verbose) vcfConcatenated = vcfConcatenated.view {
  "Variant Calling output:\n\
  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: [${it[3]}, ${it[2]}]\n\
  File  : ${it[4].fileName}"
}

process RunStrelka {
  tag {idSampleTumor + "_vs_" + idSampleNormal}

  publishDir directoryMap.strelka, mode: 'copy'

  input:
    set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from bamsForStrelka
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set val("strelka"), idPatient, idSampleNormal, idSampleTumor, file("*.vcf.gz"), file("*.vcf.gz.tbi") into strelkaOutput

  when: 'strelka' in tools

  script:
  """
  \$STRELKA_INSTALL_PATH/bin/configureStrelkaSomaticWorkflow.py \
  --tumor $bamTumor \
  --normal $bamNormal \
  --referenceFasta $genomeFile \
  --runDir Strelka

  python Strelka/runWorkflow.py -m local -j $task.cpus

  mv Strelka/results/variants/somatic.indels.vcf.gz Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_indels.vcf.gz
  mv Strelka/results/variants/somatic.indels.vcf.gz.tbi Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_indels.vcf.gz.tbi
  mv Strelka/results/variants/somatic.snvs.vcf.gz Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_snvs.vcf.gz
  mv Strelka/results/variants/somatic.snvs.vcf.gz.tbi Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_snvs.vcf.gz.tbi
  """
}

if (verbose) strelkaOutput = strelkaOutput.view {
  "Variant Calling output:\n\
  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: [${it[3]}, ${it[2]}]\n\
  Files : ${it[4].fileName}\n\
  Index : ${it[5].fileName}"
}

process RunSingleStrelka {
  tag {idSample}

  publishDir directoryMap.strelka, mode: 'copy'

  input:
    set idPatient, status, idSample, file(bam), file(bai) from bamsForSingleStrelka
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex
    ])

  output:
    set val("singlestrelka"), idPatient, idSample,  file("*.vcf.gz"), file("*.vcf.gz.tbi") into singleStrelkaOutput

  when: 'strelka' in tools

  script:
  """
  \$STRELKA_INSTALL_PATH/bin/configureStrelkaGermlineWorkflow.py \
  --bam $bam \
  --referenceFasta $genomeFile \
  --runDir Strelka

  python Strelka/runWorkflow.py -m local -j $task.cpus

  mv Strelka/results/variants/genome.*.vcf.gz Strelka_${idSample}_genome.vcf.gz
  mv Strelka/results/variants/genome.*.vcf.gz.tbi Strelka_${idSample}_genome.vcf.gz.tbi
  mv Strelka/results/variants/variants.vcf.gz Strelka_${idSample}_variants.vcf.gz
  mv Strelka/results/variants/variants.vcf.gz.tbi Strelka_${idSample}_variants.vcf.gz.tbi
  """
}

if (verbose) singleStrelkaOutput = singleStrelkaOutput.view {
  "Variant Calling output:\n\
  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: ${it[2]}\n\
  Files : ${it[3].fileName}\n\
  Index : ${it[4].fileName}"
}

process RunManta {
  tag {idSampleTumor + "_vs_" + idSampleNormal}

  publishDir directoryMap.manta, mode: 'copy'

  input:
    set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from bamsForManta
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex
    ])

  output:
    set val("manta"), idPatient, idSampleNormal, idSampleTumor, file("*.vcf.gz"), file("*.vcf.gz.tbi") into mantaOutput

  when: 'manta' in tools

  script:
  """
  \$MANTA_INSTALL_PATH/bin/configManta.py \
  --normalBam $bamNormal \
  --tumorBam $bamTumor \
  --reference $genomeFile \
  --runDir Manta

  python Manta/runWorkflow.py -m local -j $task.cpus

  mv Manta/results/variants/candidateSmallIndels.vcf.gz Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSmallIndels.vcf.gz
  mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSmallIndels.vcf.gz.tbi
  mv Manta/results/variants/candidateSV.vcf.gz Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSV.vcf.gz
  mv Manta/results/variants/candidateSV.vcf.gz.tbi Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSV.vcf.gz.tbi
  mv Manta/results/variants/diploidSV.vcf.gz Manta_${idSampleTumor}_vs_${idSampleNormal}.diploidSV.vcf.gz
  mv Manta/results/variants/diploidSV.vcf.gz.tbi Manta_${idSampleTumor}_vs_${idSampleNormal}.diploidSV.vcf.gz.tbi
  mv Manta/results/variants/somaticSV.vcf.gz Manta_${idSampleTumor}_vs_${idSampleNormal}.somaticSV.vcf.gz
  mv Manta/results/variants/somaticSV.vcf.gz.tbi Manta_${idSampleTumor}_vs_${idSampleNormal}.somaticSV.vcf.gz.tbi
  """
}

if (verbose) mantaOutput = mantaOutput.view {
  "Variant Calling output:\n\
  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: [${it[3]}, ${it[2]}]\n\
  Files : ${it[4].fileName}\n\
  Index : ${it[5].fileName}"
}

process RunSingleManta {
  tag {status == 0 ? idSample + " - Single Diploid" : idSample + " - Tumor-Only"}

  publishDir directoryMap.manta, mode: 'copy'

  input:
    set idPatient, status, idSample, file(bam), file(bai) from bamsForSingleManta
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex
    ])

  output:
    set val("singlemanta"), idPatient, idSample,  file("*.vcf.gz"), file("*.vcf.gz.tbi") into singleMantaOutput

  when: 'manta' in tools

  script:
  if ( status == 0 ) // If Normal Sample
  """
  \$MANTA_INSTALL_PATH/bin/configManta.py \
  --bam $bam \
  --reference $genomeFile \
  --runDir Manta

  python Manta/runWorkflow.py -m local -j $task.cpus

  mv Manta/results/variants/candidateSmallIndels.vcf.gz Manta_${idSample}.candidateSmallIndels.vcf.gz
  mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi Manta_${idSample}.candidateSmallIndels.vcf.gz.tbi
  mv Manta/results/variants/candidateSV.vcf.gz Manta_${idSample}.candidateSV.vcf.gz
  mv Manta/results/variants/candidateSV.vcf.gz.tbi Manta_${idSample}.candidateSV.vcf.gz.tbi
  mv Manta/results/variants/diploidSV.vcf.gz Manta_${idSample}.diploidSV.vcf.gz
  mv Manta/results/variants/diploidSV.vcf.gz.tbi Manta_${idSample}.diploidSV.vcf.gz.tbi
  """
  else  // Tumor Sample
  """
  \$MANTA_INSTALL_PATH/bin/configManta.py \
  --tumorBam $bam \
  --reference $genomeFile \
  --runDir Manta

  python Manta/runWorkflow.py -m local -j $task.cpus

  mv Manta/results/variants/candidateSmallIndels.vcf.gz Manta_${idSample}.candidateSmallIndels.vcf.gz
  mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi Manta_${idSample}.candidateSmallIndels.vcf.gz.tbi
  mv Manta/results/variants/candidateSV.vcf.gz Manta_${idSample}.candidateSV.vcf.gz
  mv Manta/results/variants/candidateSV.vcf.gz.tbi Manta_${idSample}.candidateSV.vcf.gz.tbi
  mv Manta/results/variants/tumorSV.vcf.gz Manta_${idSample}.tumorSV.vcf.gz
  mv Manta/results/variants/tumorSV.vcf.gz.tbi Manta_${idSample}.tumorSV.vcf.gz.tbi
  """
}

if (verbose) singleMantaOutput = singleMantaOutput.view {
  "Variant Calling output:\n\
  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: ${it[2]}\n\
  Files : ${it[3].fileName}\n\
  Index : ${it[4].fileName}"
}

// Run commands and code from Malin Larsson
// Based on Jesper Eisfeldt's code
process RunAlleleCount {
  tag {idSample}

  input:
    set idPatient, status, idSample, file(bam), file(bai) from bamsForAscat
    set file(acLoci), file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.acLoci,
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set idPatient, status, idSample, file("${idSample}.alleleCount") into alleleCountOutput

  when: 'ascat' in tools

  script:
  """
  alleleCounter -l $acLoci -r $genomeFile -b $bam -o ${idSample}.alleleCount;
  """
}

alleleCountNormal = Channel.create()
alleleCountTumor = Channel.create()

alleleCountOutput
  .choice(alleleCountTumor, alleleCountNormal) {it[1] == 0 ? 1 : 0}

alleleCountOutput = alleleCountNormal.combine(alleleCountTumor)

alleleCountOutput = alleleCountOutput.map {
  idPatientNormal, statusNormal, idSampleNormal, alleleCountNormal,
  idPatientTumor,  statusTumor,  idSampleTumor,  alleleCountTumor ->
  [idPatientNormal, idSampleNormal, idSampleTumor, alleleCountNormal, alleleCountTumor]
}

// R script from Malin Larssons bitbucket repo:
// https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
process RunConvertAlleleCounts {
  tag {idSampleTumor + "_vs_" + idSampleNormal}

  publishDir directoryMap.ascat, mode: 'copy'

  input:
    set idPatient, idSampleNormal, idSampleTumor, file(alleleCountNormal), file(alleleCountTumor) from alleleCountOutput

  output:
    set idPatient, idSampleNormal, idSampleTumor, file("${idSampleNormal}.BAF"), file("${idSampleNormal}.LogR"), file("${idSampleTumor}.BAF"), file("${idSampleTumor}.LogR") into convertAlleleCountsOutput

  when: 'ascat' in tools

  script:
  gender = patientGenders[idPatient]
  """
  convertAlleleCounts.r $idSampleTumor $alleleCountTumor $idSampleNormal $alleleCountNormal $gender
  """
}

// R scripts from Malin Larssons bitbucket repo:
// https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
process RunAscat {
  tag {idSampleTumor + "_vs_" + idSampleNormal}

  publishDir directoryMap.ascat, mode: 'copy'

  input:
    set idPatient, idSampleNormal, idSampleTumor, file(bafNormal), file(logrNormal), file(bafTumor), file(logrTumor) from convertAlleleCountsOutput

  output:
    set val("ascat"), idPatient, idSampleNormal, idSampleTumor, file("${idSampleTumor}.*.{png,txt}") into ascatOutput

  when: 'ascat' in tools

  script:
  """
  run_ascat.r $bafTumor $logrTumor $bafNormal $logrNormal $idSampleTumor $baseDir
  """
}

if (verbose) ascatOutput = ascatOutput.view {
  "Variant Calling output:\n\
  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: [${it[3]}, ${it[2]}]\n\
  Files : [${it[4].fileName}]"
}

vcfToAnnotate = Channel.create()
vcfNotToAnnotate = Channel.create()

if (step == 'annotate' && annotateVCF == []) {
  Channel.empty().mix(
    Channel.fromPath('VariantCalling/HaplotypeCaller/*.vcf.gz')
      .flatten().map{vcf -> ['haplotypecaller',vcf]},
    Channel.fromPath('VariantCalling/Manta/*SV.vcf.gz')
      .flatten().map{vcf -> ['manta',vcf]},
    Channel.fromPath('VariantCalling/MuTect1/*.vcf.gz')
      .flatten().map{vcf -> ['mutect1',vcf]},
    Channel.fromPath('VariantCalling/MuTect2/*.vcf.gz')
      .flatten().map{vcf -> ['mutect2',vcf]},
    Channel.fromPath('VariantCalling/Strelka/*{somatic,variants}*.vcf.gz')
      .flatten().map{vcf -> ['strelka',vcf]}
  ).choice(vcfToAnnotate, vcfNotToAnnotate) { annotateTools == [] || (annotateTools != [] && it[0] in annotateTools) ? 0 : 1 }

} else if (step == 'annotate' && annotateTools == [] && annotateVCF != []) {
  list = ""
  annotateVCF.each{ list += ",$it" }
  list = list.substring(1)
  if (StringUtils.countMatches("$list", ",") == 0) vcfToAnnotate = Channel.fromPath("$list")
    .map{vcf -> ['userspecified',vcf]}
  else vcfToAnnotate = Channel.fromPath("{$list}")
    .map{vcf -> ['userspecified',vcf]}

} else if (step != 'annotate') {
  vcfConcatenated
    .choice(vcfToAnnotate, vcfNotToAnnotate) { it[0] == 'gvcf-hc' || it[0] == 'freebayes' ? 1 : 0 }

  (strelkaIndels, strelkaSNVS) = strelkaOutput.into(2)
  (mantaSomaticSV, mantaDiploidSV) = mantaOutput.into(2)
  vcfToAnnotate = vcfToAnnotate.map {
    variantcaller, idPatient, idSampleNormal, idSampleTumor, vcf ->
    [variantcaller, vcf]
  }.mix(
    mantaDiploidSV.map {
      variantcaller, idPatient, idSampleNormal, idSampleTumor, vcf, tbi ->
      [variantcaller, vcf[2]]
    },
    mantaSomaticSV.map {
      variantcaller, idPatient, idSampleNormal, idSampleTumor, vcf, tbi ->
      [variantcaller, vcf[3]]
    },
    singleStrelkaOutput.map {
      variantcaller, idPatient, idSample, vcf, tbi ->
      [variantcaller, vcf[1]]
    },
    singleMantaOutput.map {
      variantcaller, idPatient, idSample, vcf, tbi ->
      [variantcaller, vcf[2]]
    },
    strelkaIndels.map {
      variantcaller, idPatient, idSampleNormal, idSampleTumor, vcf, tbi ->
      [variantcaller, vcf[0]]
    },
    strelkaSNVS.map {
      variantcaller, idPatient, idSampleNormal, idSampleTumor, vcf, tbi ->
      [variantcaller, vcf[1]]
    })
} else exit 1, "specify only tools or files to annotate, bot both"

vcfNotToAnnotate.close()

(vcfForBCFtools, vcfForSnpeff, vcfForVep) = vcfToAnnotate.into(3)

process RunBcftoolsStats {
  tag {vcf}

  publishDir directoryMap.bcftoolsStats, mode: 'copy'

  input:
    set variantCaller, file(vcf) from vcfForBCFtools

  output:
    file ("${vcf.baseName}.bcf.tools.stats.out") into bcfReport

  when: reports

  script:
  """
  bcftools stats $vcf > ${vcf.baseName}.bcf.tools.stats.out
  """
}

if (verbose) bcfReport = bcfReport.view {
  "BCFTools stats report:\n\
  File  : [${it.fileName}]"
}

process RunSnpeff {
  tag {vcf}

  publishDir directoryMap.snpeff, mode: 'copy'

  input:
    set variantCaller, file(vcf) from vcfForSnpeff
    val snpeffDb from Channel.value(params.genomes[params.genome].snpeffDb)

  output:
    set file("${vcf.baseName}.snpEff.ann.vcf"), file("${vcf.baseName}.snpEff.genes.txt"), file("${vcf.baseName}.snpEff.csv"), file("${vcf.baseName}.snpEff.summary.html") into snpeffReport

  when: 'snpeff' in tools

  script:
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar \$SNPEFF_HOME/snpEff.jar \
  $snpeffDb \
  -csvStats ${vcf.baseName}.snpEff.csv \
  -nodownload \
  -cancer \
  -v \
  ${vcf} \
  > ${vcf.baseName}.snpEff.ann.vcf

  mv snpEff_summary.html ${vcf.baseName}.snpEff.summary.html
  """
}

if (verbose) snpeffReport = snpeffReport.view {
  "snpEff report:\n\
  File  : ${it.fileName}"
}

process RunVEP {
  tag {vcf}

  publishDir directoryMap.vep, mode: 'copy'

  input:
    set variantCaller, file(vcf) from vcfForVep

  output:
    set file("${vcf.baseName}.vep.ann.vcf"), file("${vcf.baseName}.vep.summary.html") into vepReport

  when: 'vep' in tools

  script:
  genome = params.genome == 'smallGRCh37' ? 'GRCh37' : params.genome
  """
  vep \
  -i $vcf \
  -o ${vcf.baseName}.vep.ann.vcf \
  --stats_file ${vcf.baseName}.vep.summary.html \
  --cache \
  --everything \
  --filter_common \
  --format vcf \
  --offline \
  --pick \
  --total_length \
  --vcf
  """
}

if (verbose) vepReport = vepReport.view {
  "VEP report:\n\
  Files : ${it.fileName}"
}

process GenerateMultiQCconfig {
  publishDir directoryMap.multiQC, mode: 'copy'

  input:

  output:
  file("multiqc_config.yaml") into multiQCconfig

  when: reports

  script:
  annotateToolString = annotateTools ? "- Annotate on : ${annotateTools.join(", ")}" : ''
  annotateVCFstring = annotateVCF ? "- Annotate on : ${annotateVCF.join(", ")}" : ''
  tsvString = step != 'annotate' ? "- TSV file: ${tsvFile}" : ''
  """
  touch multiqc_config.yaml
  echo "custom_logo: $baseDir/doc/images/CAW_logo.png" >> multiqc_config.yaml
  echo "custom_logo_url: http://opensource.scilifelab.se/projects/caw" >> multiqc_config.yaml
  echo "custom_logo_title: 'Cancer Analysis Workflow'" >> multiqc_config.yaml
  echo "report_header_info:" >> multiqc_config.yaml
  echo "- CAW version: $version" >> multiqc_config.yaml
  echo "- Contact Name: ${params.callName}" >> multiqc_config.yaml
  echo "- Contact E-mail: ${params.contactMail}" >> multiqc_config.yaml
  echo "- Command Line: ${workflow.commandLine}" >> multiqc_config.yaml
  echo "- Directory: ${workflow.launchDir}" >> multiqc_config.yaml
  echo ${tsvString} >> multiqc_config.yaml
    echo "- Genome: "${params.genome} >> multiqc_config.yaml
  echo "- Step: "${step} >> multiqc_config.yaml
  echo "- Tools: "${tools.join(", ")} >> multiqc_config.yaml
  echo ${annotateToolString} >> multiqc_config.yaml
  echo ${annotateVCFstring} >> multiqc_config.yaml
  echo "  acLoci      : $referenceMap.acLoci" >> multiqc_config.yaml
  echo "  bwaIndex    : "${referenceMap.bwaIndex.join(", ")} >> multiqc_config.yaml
  echo "  cosmic      : $referenceMap.cosmic" >> multiqc_config.yaml
  echo "  cosmicIndex : $referenceMap.cosmicIndex" >> multiqc_config.yaml
  echo "  dbsnp       : $referenceMap.dbsnp" >> multiqc_config.yaml
  echo "  dbsnpIndex  : $referenceMap.dbsnpIndex" >> multiqc_config.yaml
  echo "  genomeDict  : $referenceMap.genomeDict" >> multiqc_config.yaml
  echo "  genomeFile  : $referenceMap.genomeFile" >> multiqc_config.yaml
  echo "  genomeIndex : $referenceMap.genomeIndex" >> multiqc_config.yaml
  echo "  intervals   : $referenceMap.intervals" >> multiqc_config.yaml
  echo "  knownIndels : "${referenceMap.knownIndels.join(", ")} >> multiqc_config.yaml
  echo "  knownIndelsIndex: "${referenceMap.knownIndelsIndex.join(", ")} >> multiqc_config.yaml
  echo "  snpeffDb    : ${params.genomes[params.genome].snpeffDb}" >> multiqc_config.yaml
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
    Channel.fromPath('Reports/{BCFToolsStats,MarkDuplicates,SamToolsStats}/*'),
    Channel.fromPath('Reports/{bamQC,FastQC}/*/*'),
    bamQCreport,
    bcfReport,
    fastQCreport,
    markDuplicatesReport,
    multiQCconfig,
    samtoolsStatsReport,
    snpeffReport,
    vepReport
  ).collect()

process RunMultiQC {
  publishDir directoryMap.multiQC, mode: 'copy'

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

def cawMessage() {
  // Display CAW message
  log.info "CANCER ANALYSIS WORKFLOW ~ $version - " + this.grabRevision() + (workflow.commitId ? " [$workflow.commitId]" : "")
}

def checkFileExtension(it, extension) {
  // Check file extension
  if (!it.toString().toLowerCase().endsWith(extension.toLowerCase())) {
    exit 1, "File: $it has the wrong extension: $extension see --help for more information"
  }
}

def checkParameterExistence(it, list) {
  // Check parameter existence
  if (!list.contains(it)) {
    println("Unknown parameter: $it")
    return false
  }
  return true
}

def checkParameterList(list, realList) {
  // Loop through all parameters to check their existence and spelling
  return list.every{ checkParameterExistence(it, realList) }
}

def checkParamReturnFile(item) {
  params."$item" = params.genomes[params.genome]."$item"
  return file(params."$item")
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
    'no-GVCF',
    'no-reports',
    'noGVCF',
    'noReports',
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

def checkReferenceMap(referenceMap) {
  // Loop through all the references files to check their existence
  referenceMap.every {
    referenceFile, fileToCheck ->
    checkRefExistence(referenceFile, fileToCheck)
  }
}

def checkRefExistence(referenceFile, fileToCheck) {
  if (fileToCheck instanceof List) {
    return fileToCheck.every{ checkRefExistence(referenceFile, it) }
  }
  def f = file(fileToCheck)
  if (f instanceof List && f.size() > 0) {
    // this is an expanded wildcard: we can assume all files exist
    return true
  } else if (!f.exists()) {
    log.info  "Missing references: $referenceFile $fileToCheck"
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
    'nonRealigned'     : 'Preprocessing/NonRealigned',
    'nonRecalibrated'  : 'Preprocessing/NonRecalibrated',
    'recalibrated'     : 'Preprocessing/Recalibrated',
    'bamQC'            : 'Reports/bamQC',
    'bcftoolsStats'    : 'Reports/BCFToolsStats',
    'fastQC'           : 'Reports/FastQC',
    'markDuplicatesQC' : 'Reports/MarkDuplicates',
    'multiQC'          : 'Reports/MultiQC',
    'samtoolsStats'    : 'Reports/SamToolsStats',
    'ascat'            : 'VariantCalling/Ascat',
    'freebayes'        : 'VariantCalling/FreeBayes',
    'haplotypecaller'  : 'VariantCalling/HaplotypeCaller',
    'gvcf-hc'          : 'VariantCalling/HaplotypeCallerGVCF',
    'manta'            : 'VariantCalling/Manta',
    'mutect1'          : 'VariantCalling/MuTect1',
    'mutect2'          : 'VariantCalling/MuTect2',
    'strelka'          : 'VariantCalling/Strelka',
    'snpeff'           : 'Annotation/SnpEff',
    'vep'              : 'Annotation/VEP'
  ]
}

def defineReferenceMap() {
  if (!(params.genome in params.genomes)) {
    exit 1, "Genome $params.genome not found in configuration"
  }
  return [
    // loci file for ascat
    'acLoci'           : checkParamReturnFile("acLoci"),
    'dbsnp'            : checkParamReturnFile("dbsnp"),
    'dbsnpIndex'       : checkParamReturnFile("dbsnpIndex"),
    // cosmic VCF with VCF4.1 header
    'cosmic'           : checkParamReturnFile("cosmic"),
    'cosmicIndex'      : checkParamReturnFile("cosmicIndex"),
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
    'annotate',
    'mapping',
    'realign',
    'recalibrate',
    'variantcalling'
  ]
}

def defineToolList() {
  return [
    'ascat',
    'freebayes',
    'haplotypecaller',
    'manta',
    'mutect1',
    'mutect2',
    'snpeff',
    'strelka',
    'vep'
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
      def fastqFile1 = workflow.commitId && params.test ? returnFile("$workflow.projectDir/${list[5]}") : returnFile("${list[5]}")
      def fastqFile2 = workflow.commitId && params.test ? returnFile("$workflow.projectDir/${list[6]}") : returnFile("${list[6]}")

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
    .ifEmpty { error "No directories found matching pattern '$pattern'" }
    .subscribe onNext: { sampleDir ->
      // the last name of the sampleDir is assumed to be a unique sample id
      sampleId = sampleDir.getFileName().toString()

      for (path1 in file("${sampleDir}/**_R1_*.fastq.gz")) {
        assert path1.getName().contains('_R1_')
        path2 = file(path1.toString().replace('_R1_', '_R2_'))
        if (!path2.exists()) {
            error "Path '${path2}' not found"
        }
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

def generateIntervalsForVC(bams, intervals) {

  def (bamsNew, bamsForVC) = bams.into(2)
  def (intervalsNew, vcIntervals) = intervals.into(2)
  def bamsForVCNew = bamsForVC.combine(vcIntervals)
  return [bamsForVCNew, bamsNew, intervalsNew]
}

def grabRevision() {
  // Return the same string executed from github or not
  return workflow.revision ?: workflow.commitId ?: workflow.scriptId.substring(0,10)
}

def helpMessage() {
  // Display help message
  this.cawMessage()
  log.info "    Usage:"
  log.info "       nextflow run SciLifeLab/CAW --sample <file.tsv> [--step STEP] [--tools TOOL[,TOOL]] --genome <Genome>"
  log.info "       nextflow run SciLifeLab/CAW --sampleDir <Directory> [--step STEP] [--tools TOOL[,TOOL]] --genome <Genome>"
  log.info "       nextflow run SciLifeLab/CAW --test [--step STEP] [--tools TOOL[,TOOL]] --genome <Genome>"
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
  log.info "         variantcalling (will start workflow with recalibrated BAM files)"
  log.info "         annotate (will annotate Variant Calling output."
  log.info "         By default it will try to annotate all available vcfs."
  log.info "         Use with --annotateTools or --annotateVCF to specify what to annotate"
  log.info "    --noReports"
  log.info "       Disable QC tools and MultiQC to generate a HTML report"
  log.info "    --tools"
  log.info "       Option to configure which tools to use in the workflow."
  log.info "         Different tools to be separated by commas."
  log.info "       Possible values are:"
  log.info "         mutect1 (use MuTect1 for VC)"
  log.info "         mutect2 (use MuTect2 for VC)"
  log.info "         freebayes (use FreeBayes for VC)"
  log.info "         strelka (use Strelka for VC)"
  log.info "         haplotypecaller (use HaplotypeCaller for normal bams VC)"
  log.info "         manta (use Manta for SV)"
  log.info "         ascat (use Ascat for CNV)"
  log.info "         snpeff (use snpEff for Annotation of Variants)"
  log.info "         vep (use VEP for Annotation of Variants)"
  log.info "    --annotateTools"
  log.info "       Option to configure which tools to annotate."
  log.info "         Different tools to be separated by commas."
  log.info "       Possible values are:"
  log.info "         haplotypecaller (Annotate HaplotypeCaller output)"
  log.info "         manta (Annotate Manta output)"
  log.info "         mutect1 (Annotate MuTect1 output)"
  log.info "         mutect2 (Annotate MuTect2 output)"
  log.info "         strelka (Annotate Strelka output)"
  log.info "    --annotateVCF"
  log.info "       Option to configure which vcf to annotate."
  log.info "         Different vcf to be separated by commas."
  log.info "    --genome <Genome>"
  log.info "       Use a specific genome version."
  log.info "       Possible values are:"
  log.info "         GRCh37"
  log.info "         GRCh38 (Default)"
  log.info "         smallGRCh37 (Use a small reference (Tests only))"
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
  log.info "Command Line: $workflow.commandLine"
  log.info "Project Dir : $workflow.projectDir"
  log.info "Launch Dir  : $workflow.launchDir"
  log.info "Work Dir    : $workflow.workDir"
  if (step != 'annotate') log.info "TSV file    : $tsvFile"
  log.info "Genome      : " + params.genome
  log.info "Step        : " + step
  if (tools) log.info "Tools       : " + tools.join(', ')
  if (annotateTools) log.info "Annotate on : " + annotateTools.join(', ')
  if (annotateVCF) log.info "VCF files   : " +annotateVCF.join(',\n    ')
  log.info "Reference files used:"
  log.info "  acLoci      : $referenceMap.acLoci"
  log.info "  bwaIndex    : " + referenceMap.bwaIndex.join(',\n    ')
  log.info "  cosmic      : $referenceMap.cosmic"
  log.info "  cosmicIndex : $referenceMap.cosmicIndex"
  log.info "  dbsnp       : $referenceMap.dbsnp"
  log.info "  dbsnpIndex  : $referenceMap.dbsnpIndex"
  log.info "  genomeDict  : $referenceMap.genomeDict"
  log.info "  genomeFile  : $referenceMap.genomeFile"
  log.info "  genomeIndex : $referenceMap.genomeIndex"
  log.info "  intervals   : $referenceMap.intervals"
  log.info "  knownIndels : " + referenceMap.knownIndels.join(',\n    ')
  log.info "  knownIndelsIndex: " + referenceMap.knownIndelsIndex.join(',\n    ')
  log.info "  snpeffDb    : ${params.genomes[params.genome].snpeffDb}"
}

def nextflowMessage() {
  // Nextflow message (version + build)
  log.info "N E X T F L O W  ~  version $workflow.nextflow.version $workflow.nextflow.build"
}

def returnFile(it) {
  // return file if it exists
  final f = file(it)
  if (!f.exists()) {
    exit 1, "Missing file in TSV file: $it, see --help for more information"
  }
  return f
}

def returnStatus(it) {
  // Return status if it's correct
  // Status should be only 0 or 1
  // 0 being normal
  // 1 being tumor (or relapse or anything that is not normal...)
  if (!(it in [0, 1])) {
    exit 1, "Status is not recognized in TSV file: $it, see --help for more information"
  }
  return it
}

def returnTSV(it, number) {
  // return TSV if it has the correct number of items in row
  if (it.size() != number) {
    exit 1, "Malformed row in TSV file: $it, see --help for more information"
  }
  return it
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
