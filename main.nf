#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
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
 - MapReads - Map reads
 - MergeBams - Merge BAMs if multilane samples
 - MarkDuplicates - Mark Duplicates
 - CreateIntervals - Create Intervals
 - RealignBams - Realign Bams as T/N pair
 - CreateRecalibrationTable - Create Recalibration Table
 - RecalibrateBam - Recalibrate Bam
 - RunSamtoolsStats - Run Samtools stats on recalibrated BAM files
 - RunHaplotypecaller - Run HaplotypeCaller for GermLine Variant Calling (Parrallelized processes)
 - RunMutect1 - Run MuTect1 for Variant Calling (Parrallelized processes)
 - RunMutect2 - Run MuTect2 for Variant Calling (Parrallelized processes)
 - RunFreeBayes - Run FreeBayes for Variant Calling (Parrallelized processes)
 - RunVardict - Run VarDict for Variant Calling (Parrallelized processes)
 - ConcatVCF - Merge results from HaplotypeCaller, MuTect1, MuTect2 and VarDict
 - RunStrelka - Run Strelka for Variant Calling
 - RunManta - Run Manta for Structural Variant Calling
 - RunAlleleCount - Run AlleleCount to prepare for ASCAT
 - RunConvertAlleleCounts - Run convertAlleleCounts to prepare for ASCAT
 - RunAscat - Run ASCAT for CNV
 - RunSnpeff - Run snpEff for annotation of vcf files
 - RunMultiQC - Run MultiQC for report and QC
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

revision = grabGitRevision()
testFile = ''
testSteps = []
version = '1.0'
workflowSteps = []

if (!checkUppmaxProject()) {exit 1, 'No UPPMAX project ID found! Use --project <UPPMAX Project ID>'}

if (params.help) {
  help_message(version, revision)
  exit 1
}
if (params.version) {
  version_message(version, revision)
  exit 1
}

workflowSteps = params.steps ? params.steps.split(',').collect {it.trim()} : ''

directoryMap = defineDirectoryMap()
referenceMap = defineReferenceMap()
stepList = defineStepList()
verbose = params.verbose ? true : false

if (!checkReferenceMap(referenceMap)) {exit 1, 'Missing Reference file(s), see --help for more information'}
if (!checkStepList(workflowSteps,stepList)) {exit 1, 'Unknown step(s), see --help for more information'}

if (params.test) {
  test = true
  referenceMap.put("intervals", "$workflow.projectDir/repeats/tiny.list")
  testFile = 'preprocessing' in workflowSteps ? file("$workflow.projectDir/data/tsv/tiny.tsv") : Channel.empty()
  testFile = 'realign' in workflowSteps ? file("$workflow.launchDir/${directoryMap['nonRealigned']}/nonRealigned.tsv") : testFile
  testFile = 'recalibrate' in workflowSteps ? file("$workflow.launchDir/${directoryMap['nonRecalibrated']}/nonRecalibrated.tsv") : testFile
  testFile = 'skipPreprocessing' in workflowSteps ? file("$workflow.launchDir/${directoryMap['recalibrated']}/recalibrated.tsv") : testFile
} else {test = false}

if (!checkSteps(workflowSteps)) {exit 1, 'Please choose only one step between preprocessing, realign, recalibrate and skipPreprocessing, see --help for more information'}

// Extract and verify content of TSV file

if ((!params.sample) && !(test)) {exit 1, 'Missing TSV file, see --help for more information'}

tsvFile = test ? testFile : file(params.sample)

fastqFiles = 'preprocessing' in workflowSteps ? extractFastqFiles(tsvFile) : Channel.empty()
bamFiles = 'realign' in workflowSteps ? extractBamFiles(tsvFile) : Channel.empty()
bamFiles = 'recalibrate' in workflowSteps ? extractRecalibrationTables(tsvFile) : bamFiles
bamFiles = 'skipPreprocessing' in workflowSteps ? extractBamFiles(tsvFile) : bamFiles

verbose ? fastqFiles = fastqFiles.view {"FASTQ files to preprocess: $it"} : ''
verbose ? bamFiles = bamFiles.view {"BAM files to process: $it"} : ''
start_message(version, revision)

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

(fastqFiles, fastqFilesforFastQC) = fastqFiles.into(2)

verbose ? fastqFilesforFastQC = fastqFilesforFastQC.view {"FASTQ files for FastQC: $it"} : ''

process RunFastQC {
  tag {idPatient + "-" + idRun}

  publishDir directoryMap['FastQC'], mode: 'copy'

  input:
    set idPatient, gender, status, idSample, idRun, file(fastqFile1), file(fastqFile2) from fastqFilesforFastQC

  output:
    file "*_fastqc.{zip,html}" into fastQCreport

  when: 'preprocessing' in workflowSteps && 'MultiQC' in workflowSteps

  script:
  """
  fastqc -q $fastqFile1 $fastqFile2
  """
}

verbose ? fastQCreport = fastQCreport.view {"FastQC report: $it"} : ''

process MapReads {
  tag {idPatient + "-" + idRun}

  input:
    set idPatient, gender, status, idSample, idRun, file(fastqFile1), file(fastqFile2) from fastqFiles
    file genomeAmb from file(referenceMap['genomeAmb'])
    file genomeAnn from file(referenceMap['genomeAnn'])
    file genomeBwt from file(referenceMap['genomeBwt'])
    file genomeFile from file(referenceMap['genomeFile'])
    file genomeIndex from file(referenceMap['genomeIndex'])
    file genomePac from file(referenceMap['genomePac'])
    file genomeSa from file(referenceMap['genomeSa'])

  output:
    set idPatient, gender, status, idSample, idRun, file("${idRun}.bam") into bam

  when: 'preprocessing' in workflowSteps

  script:
  readGroup="@RG\\tID:$idRun\\tSM:$idSample\\tLB:$idSample\\tPL:illumina"
  """
  set -eo pipefail
  bwa mem -R \"$readGroup\" -B 3 -t $task.cpus -M \
  $genomeFile $fastqFile1 $fastqFile2 | \
  samtools sort --threads $task.cpus - > ${idRun}.bam
  """
}

verbose ? bam = bam.view {"BAM file to sort into group or single: $it"} : ''

// Sort bam whether they are standalone or should be merged
// Borrowed code from https://github.com/guigolab/chip-nf

singleBam = Channel.create()
groupedBam = Channel.create()
bam.groupTuple(by:[0,1,2,3])
  .choice(singleBam, groupedBam) {it[4].size() > 1 ? 1 : 0}
singleBam = singleBam.map {
  idPatient, gender, status, idSample, idRun, bam ->
  [idPatient, gender, status, idSample, bam]
}

verbose ? groupedBam = groupedBam.view {"Grouped BAMs to merge: $it"} : ''

process MergeBams {
  tag {idPatient + "-" + idSample}

  input:
    set idPatient, gender, status, idSample, idRun, file(bam) from groupedBam

  output:
    set idPatient, gender, status, idSample, file("${idSample}.bam") into mergedBam

  when: 'preprocessing' in workflowSteps

  script:
  """
  samtools merge --threads $task.cpus ${idSample}.bam $bam
  """
}

verbose ? singleBam = singleBam.view {"Single BAM: $it"} : ''
verbose ? mergedBam = mergedBam.view {"Merged BAM: $it"} : ''
mergedBam = mergedBam.mix(singleBam)
verbose ? mergedBam = mergedBam.view {"BAM for MarkDuplicates: $it"} : ''

process MarkDuplicates {
  tag {idPatient + "-" + idSample}

  publishDir '.', saveAs: { it == "${bam}.metrics" ? "${directoryMap['MarkDuplicatesQC']}/$it" : "${directoryMap['nonRealigned']}/$it" }, mode: 'copy'

  input:
    set idPatient, gender, status, idSample, file(bam) from mergedBam

  output:
    set idPatient, gender, val("${idSample}_${status}"), file("${idSample}_${status}.md.bam"), file("${idSample}_${status}.md.bai") into duplicates
    set idPatient, gender, status, idSample, val("${idSample}_${status}.md.bam"), val("${idSample}_${status}.md.bai") into markDuplicatesTSV
    file ("${bam}.metrics") into markDuplicatesReport

  when: 'preprocessing' in workflowSteps

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
markDuplicatesTSV.map { idPatient, gender, status, idSample, bam, bai ->
  "$idPatient\t$gender\t$status\t$idSample\t${directoryMap['nonRealigned']}/$bam\t${directoryMap['nonRealigned']}/$bai\n"
}.collectFile( name: 'nonRealigned.tsv', sort: true, storeDir: directoryMap['nonRealigned'])

// Create intervals for realignement using both tumor+normal as input
// Group the marked duplicates BAMs for intervals and realign by idPatient
// Grouping also by gender, to make a nicer channel
duplicatesGrouped = 'preprocessing' in workflowSteps ? duplicates.groupTuple(by:[0,1]) : Channel.empty()

duplicatesGrouped = 'realign' in workflowSteps ? bamFiles.map{
  idPatient, gender, status, idSample, bam, bai ->
  [idPatient, gender, "${idSample}_${status}", bam, bai]
}.groupTuple(by:[0,1]) : duplicatesGrouped

// The duplicatesGrouped channel is duplicated
// one copy goes to the CreateIntervals process
// and the other to the RealignBams process
(duplicatesInterval, duplicatesRealign) = duplicatesGrouped.into(2)

verbose ? duplicatesInterval = duplicatesInterval.view {"BAMs for CreateIntervals: $it"} : ''
verbose ? duplicatesRealign = duplicatesRealign.view {"BAMs to phase: $it"} : ''
verbose ? markDuplicatesReport = markDuplicatesReport.view {"MarkDuplicates report: $it"} : ''

// VCF indexes are added so they will be linked, and not re-created on the fly
//  -L "1:131941-141339" \
process CreateIntervals {
  tag {idPatient}

  input:
    set idPatient, gender, idSample_status, file(bam), file(bai) from duplicatesInterval
    file genomeFile from file(referenceMap['genomeFile'])
    file genomeIndex from file(referenceMap['genomeIndex'])
    file genomeDict from file(referenceMap['genomeDict'])
    file kgIndels from file(referenceMap['kgIndels'])
    file kgIndex from file(referenceMap['kgIndex'])
    file millsIndels from file(referenceMap['millsIndels'])
    file millsIndex from file(referenceMap['millsIndex'])

  output:
    set idPatient, gender, file("${idPatient}.intervals") into intervals

  when: 'preprocessing' in workflowSteps || 'realign' in workflowSteps

  script:
  bams = bam.collect{"-I $it"}.join(' ')
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar \$GATK_HOME/GenomeAnalysisTK.jar \
  -T RealignerTargetCreator \
  $bams \
  -R $genomeFile \
  -known $kgIndels \
  -known $millsIndels \
  -nt $task.cpus \
  -XL hs37d5 \
  -XL NC_007605 \
  -o ${idPatient}.intervals
  """
}

verbose ? intervals = intervals.view {"Intervals to phase: $it"} : ''

bamsAndIntervals = duplicatesRealign
  .phase(intervals)
  .map{duplicatesRealign, intervals ->
    tuple(
      duplicatesRealign[0],
      duplicatesRealign[1],
      duplicatesRealign[2],
      duplicatesRealign[3],
      duplicatesRealign[4],
      intervals[2]
    )}

verbose ? bamsAndIntervals = bamsAndIntervals.view {"Bams and Intervals phased for RealignBams: $it"} : ''

// use nWayOut to split into T/N pair again
process RealignBams {
  tag {idPatient}

  input:
    set idPatient, gender, idSample_status, file(bam), file(bai), file(intervals) from bamsAndIntervals
    file genomeFile from file(referenceMap['genomeFile'])
    file genomeIndex from file(referenceMap['genomeIndex'])
    file genomeDict from file(referenceMap['genomeDict'])
    file kgIndels from file(referenceMap['kgIndels'])
    file kgIndex from file(referenceMap['kgIndex'])
    file millsIndels from file(referenceMap['millsIndels'])
    file millsIndex from file(referenceMap['millsIndex'])

  output:
    set idPatient, gender, file("*.real.bam"), file("*.real.bai") into realignedBam mode flatten

  when: 'preprocessing' in workflowSteps || 'realign' in workflowSteps

  script:
  bams = bam.collect{"-I $it"}.join(' ')
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar \$GATK_HOME/GenomeAnalysisTK.jar \
  -T IndelRealigner \
  $bams \
  -R $genomeFile \
  -targetIntervals $intervals \
  -known $kgIndels \
  -known $millsIndels \
  -XL hs37d5 \
  -XL NC_007605 \
  -nWayOut '.real.bam'
  """
}

realignedBam = retreiveStatus(realignedBam)

verbose ? realignedBam = realignedBam.view {"Realigned BAM to CreateRecalibrationTable: $it"} : ''

process CreateRecalibrationTable {
  tag {idPatient + "-" + idSample}

  publishDir directoryMap['nonRecalibrated'], mode: 'copy'

  input:
    set idPatient, gender, status, idSample, file(bam), file(bai) from realignedBam
    file genomeFile from file(referenceMap['genomeFile'])
    file genomeIndex from file(referenceMap['genomeIndex'])
    file genomeDict from file(referenceMap['genomeDict'])
    file dbsnp from file(referenceMap['dbsnp'])
    file dbsnpIndex from file(referenceMap['dbsnpIndex'])
    file kgIndels from file(referenceMap['kgIndels'])
    file kgIndex from file(referenceMap['kgIndex'])
    file millsIndels from file(referenceMap['millsIndels'])
    file millsIndex from file(referenceMap['millsIndex'])

  output:
    set idPatient, gender, status, idSample, file(bam), file(bai), file("${idSample}.recal.table") into recalibrationTable
    set idPatient, gender, status, idSample, val("${idSample}_${status}.md.real.bam"), val("${idSample}_${status}.md.real.bai"), val("${idSample}.recal.table") into recalibrationTableTSV

  when: 'preprocessing' in workflowSteps || 'realign' in workflowSteps

  script:
  """
  java -Xmx${task.memory.toGiga()}g \
  -Djava.io.tmpdir="/tmp" \
  -jar \$GATK_HOME/GenomeAnalysisTK.jar \
  -T BaseRecalibrator \
  -R $genomeFile \
  -I $bam \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  -knownSites $dbsnp \
  -knownSites $kgIndels \
  -knownSites $millsIndels \
  -nct $task.cpus \
  -XL hs37d5 \
  -XL NC_007605 \
  -l INFO \
  -o ${idSample}.recal.table
  """
}

// Creating a TSV file to restart from this step
recalibrationTableTSV.map { idPatient, gender, status, idSample, bam, bai, recalTable ->
  "$idPatient\t$gender\t$status\t$idSample\t${directoryMap['nonRecalibrated']}/$bam\t${directoryMap['nonRecalibrated']}/$bai\t\t${directoryMap['nonRecalibrated']}/$recalTable\n"
}.collectFile( name: 'nonRecalibrated.tsv', sort: true, storeDir: directoryMap['nonRecalibrated'])

recalibrationTable = 'recalibrate' in workflowSteps ? bamFiles : recalibrationTable

verbose ? recalibrationTable = recalibrationTable.view {"Base recalibrated table for RecalibrateBam: $it"} : ''

process RecalibrateBam {
  tag {idPatient + "-" + idSample}

  publishDir directoryMap['recalibrated'], mode: 'copy'

  input:
    set idPatient, gender, status, idSample, file(bam), file(bai), recalibrationReport from recalibrationTable
    file genomeFile from file(referenceMap['genomeFile'])
    file genomeIndex from file(referenceMap['genomeIndex'])
    file genomeDict from file(referenceMap['genomeDict'])

  output:
    set idPatient, gender, status, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bai") into recalibratedBam
    set idPatient, gender, status, idSample, val("${idSample}.recal.bam"), val("${idSample}.recal.bai") into recalibratedBamTSV

  when: 'preprocessing' in workflowSteps || 'realign' in workflowSteps || 'recalibrate' in workflowSteps

  // TODO: ditto as at the previous BaseRecalibrator step, consider using -nct 4
  script:
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar \$GATK_HOME/GenomeAnalysisTK.jar \
  -T PrintReads \
  -R $genomeFile \
  -nct $task.cpus \
  -I $bam \
  -XL hs37d5 \
  -XL NC_007605 \
  --BQSR $recalibrationReport \
  -o ${idSample}.recal.bam
  """
}

// Creating a TSV file to restart from this step
recalibratedBamTSV.map { idPatient, gender, status, idSample, bam, bai ->
  "$idPatient\t$gender\t$status\t$idSample\t${directoryMap['recalibrated']}/$bam\t${directoryMap['recalibrated']}/$bai\n"
}.collectFile( name: 'recalibrated.tsv', sort: true, storeDir: directoryMap['recalibrated'])

recalibratedBam = 'skipPreprocessing' in workflowSteps ? bamFiles : recalibratedBam

verbose ? recalibratedBam = recalibratedBam.view {"Recalibrated Bam for variant Calling: $it"} : ''

(recalibratedBam, recalibratedBamForStats) = recalibratedBam.into(2)

process RunSamtoolsStats {
  tag {idPatient + "-" + idSample}

  publishDir directoryMap['SamToolsStats'], mode: 'copy'

  input:
    set idPatient, gender, status, idSample, file(bam), file(bai) from recalibratedBamForStats

  output:
    file ("${bam}.samtools.stats.out") into recalibratedBamReport

    when: 'MultiQC' in workflowSteps

    script:
    """
    samtools stats $bam > ${bam}.samtools.stats.out
    """
}

verbose ? recalibratedBamReport = recalibratedBamReport.view {"BAM Stats: $it"} : ''

// Here we have a recalibrated bam set, but we need to separate the bam files based on patient status.
// The sample tsv config file which is formatted like: "subject status sample lane fastq1 fastq2"
// cf fastqFiles channel, I decided just to add _status to the sample name to have less changes to do.
// And so I'm sorting the channel if the sample match _0, then it's a normal sample, otherwise tumor.
// Then spread normal over tumor to get each possibilities
// ie. normal vs tumor1, normal vs tumor2, normal vs tumor3
// then copy this channel into channels for each variant calling
// I guess it will still work even if we have multiple normal samples

// separate recalibrateBams by status
bamsNormal = Channel.create()
bamsTumor = Channel.create()
recalibratedBam
  .choice(bamsTumor, bamsNormal) {it[2] =~ /^0$/ ? 1 : 0}

// Removing status because not relevant anymore
bamsNormal = bamsNormal.map { idPatient, gender, status, idSample, bam, bai -> [idPatient, gender, idSample, bam, bai] }
verbose ? bamsNormal = bamsNormal.view {"Normal Bam for variant Calling: $it"} : ''

bamsTumor = bamsTumor.map { idPatient, gender, status, idSample, bam, bai -> [idPatient, gender, idSample, bam, bai] }
verbose ? bamsTumor = bamsTumor.view {"Tumor Bam for variant Calling: $it"} : ''

// We know that MuTect2 (and other somatic callers) are notoriously slow. To speed them up we are chopping the reference into
// smaller pieces at centromeres (see repeats/centromeres.list), do variant calling by this intervals, and re-merge the VCFs.
// Since we are on a cluster, this can parallelize the variant call processes, and push down the variant call wall clock time significanlty.

// in fact we need two channels: one for the actual genomic region, and an other for names
// without ":", as nextflow is not happy with them (will report as a failed process).
// For region 1:1-2000 the output file name will be something like 1_1-2000_Sample_name.mutect2.vcf
// from the "1:1-2000" string make ["1:1-2000","1_1-2000"]

// define intervals file by --intervals
intervals = Channel.from(file(referenceMap['intervals']).readLines())
gI = intervals.map{[it,it.replaceFirst(/\:/,'_')]}

(bamsNormalTemp, bamsNormal, gI) = generateIntervalsForVC(bamsNormal, gI)
(bamsTumorTemp, bamsTumor, gI) = generateIntervalsForVC(bamsTumor, gI)

// HaplotypeCaller
bamsFHC = bamsNormalTemp.mix(bamsTumorTemp)
verbose ? bamsFHC = bamsFHC.view {"Bams with Intervals for HaplotypeCaller: $it"} : ''
if (!'HaplotypeCaller' in workflowSteps) {bamsFHC.close()}

(bamsNormalTemp, bamsNormal) = bamsNormal.into(2)
(bamsTumorTemp, bamsTumor) = bamsTumor.into(2)

bamsNormalTemp = bamsNormalTemp.map { idPatient, gender, idSample, bam, bai -> [idPatient, gender, 0, idSample, bam, bai] }
bamsTumorTemp = bamsTumorTemp.map { idPatient, gender, idSample, bam, bai -> [idPatient, gender, 1, idSample, bam, bai] }

bamsForAscat = Channel.create()
bamsForAscat = bamsNormalTemp.mix(bamsTumorTemp)
verbose ? bamsForAscat = bamsForAscat.view {"Bams for Ascat: $it"} : ''
if (!'Ascat' in workflowSteps) {bamsForAscat.close()}

bamsAll = bamsNormal.spread(bamsTumor)
// Since idPatientNormal and idPatientTumor are the same
// It's removed from bamsAll Channel (same for genderNormal)
// /!\ It is assumed that every sample are from the same patient
bamsAll = bamsAll.map {
  idPatientNormal, genderNormal, idSampleNormal, bamNormal, baiNormal, idPatientTumor, genderTumor, idSampleTumor, bamTumor, baiTumor ->
  [idPatientNormal, genderNormal, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor]
}
verbose ? bamsAll = bamsAll.view {"Mapped Recalibrated Bam for variant Calling: $it"} : ''

// MuTect1
(bamsFMT1, bamsAll, gI) = generateIntervalsForVC(bamsAll, gI)
verbose ? bamsFMT1 = bamsFMT1.view {"Bams with Intervals for MuTect1: $it"} : ''
if (!'MuTect1' in workflowSteps) {bamsFMT1.close()}

// MuTect2
(bamsFMT2, bamsAll, gI) = generateIntervalsForVC(bamsAll, gI)
verbose ? bamsFMT2 = bamsFMT2.view {"Bams with Intervals for MuTect2: $it"} : ''
if (!'MuTect2' in workflowSteps) {bamsFMT2.close()}

// FreeBayes
(bamsFFB, bamsAll, gI) = generateIntervalsForVC(bamsAll, gI)
verbose ? bamsFFB = bamsFFB.view {"Bams with Intervals for FreeBayes: $it"} : ''
if (!'FreeBayes' in workflowSteps) {bamsFFB.close()}

// VarDict
(bamsFVD, bamsAll, gI) = generateIntervalsForVC(bamsAll, gI)
verbose ? bamsFVD = bamsFVD.view {"Bams with Intervals for VarDict: $it"} : ''
if (!'VarDict' in workflowSteps) {bamsFVD.close()}

(bamsForManta, bamsForStrelka) = bamsAll.into(2)

verbose ? bamsForManta = bamsForManta.view {"Bams for Manta: $it"} : ''
if (!'Manta' in workflowSteps) {bamsForManta.close()}

verbose ? bamsForStrelka = bamsForStrelka.view {"Bams for Strelka: $it"} : ''
if (!'Strelka' in workflowSteps) {bamsForStrelka.close()}

process RunHaplotypecaller {
  tag {idPatient + "-" + idSample + "-" + gen_int}

  input:
    set idPatient, gender, idSample, file(bam), file(bai), genInt, gen_int from bamsFHC //Are these values `ped to bamNormal already?
    file genomeFile from file(referenceMap['genomeFile'])
    file genomeIndex from file(referenceMap['genomeIndex'])
    file genomeDict from file(referenceMap['genomeDict'])
    file dbsnp from file(referenceMap['dbsnp'])
    file dbsnpIndex from file(referenceMap['dbsnpIndex'])

  output:
    set val("HaplotypeCaller"), idPatient, gender, idSample, val("${gen_int}_${idSample}"), file("${gen_int}_${idSample}.vcf") into hcVCF

  when: 'HaplotypeCaller' in workflowSteps

  // both -nt and -nct removed : it is still not recommended to use more threads for HC, use scatter-gather instead
  script:
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar \$GATK_HOME/GenomeAnalysisTK.jar \
  -T HaplotypeCaller \
  -R $genomeFile \
  --dbsnp $dbsnp \
  -I $bam \
  -L \"$genInt\" \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  -XL hs37d5 \
  -XL NC_007605 \
  -o ${gen_int}_${idSample}.vcf
  """
}

hcVCF = hcVCF.map {
  variantCaller, idPatient, gender, idSample, tag, vcfFile ->
  [variantCaller, idPatient, gender, idSample, idSample, tag, vcfFile]
}.groupTuple(by:[0,1,2,3,4])

verbose ? hcVCF = hcVCF.view {"HaplotypeCaller output: $it"} : ''

process RunMutect1 {
  tag {idPatient + "-" + idSampleTumor + "-" + gen_int}

  input:
    set idPatient, gender, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), genInt, gen_int from bamsFMT1
    file genomeFile from file(referenceMap['genomeFile'])
    file genomeIndex from file(referenceMap['genomeIndex'])
    file genomeDict from file(referenceMap['genomeDict'])
    file dbsnp from file(referenceMap['dbsnp'])
    file dbsnpIndex from file(referenceMap['dbsnpIndex'])
    file cosmic from file(referenceMap['cosmic'])
    file cosmicIndex from file(referenceMap['cosmicIndex'])

  output:
    set val("MuTect1"), idPatient, gender, idSampleNormal, idSampleTumor, val("${gen_int}_${idSampleTumor}_vs_${idSampleNormal}"), file("${gen_int}_${idSampleTumor}_vs_${idSampleNormal}.vcf") into mutect1Output

  when: 'MuTect1' in workflowSteps

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
  -L \"$genInt\" \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  -XL hs37d5 \
  -XL NC_007605 \
  --out ${gen_int}_${idSampleTumor}_vs_${idSampleNormal}.call_stats.out \
  --vcf ${gen_int}_${idSampleTumor}_vs_${idSampleNormal}.vcf
  """
}

mutect1Output = mutect1Output.groupTuple(by:[0,1,2,3,4])
verbose ? mutect1Output = mutect1Output.view {"MuTect1 output: $it"} : ''

process RunMutect2 {
  tag {idPatient + "-" + idSampleTumor + "-" + gen_int}

  input:
    set idPatient, gender, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), genInt, gen_int from bamsFMT2
    file genomeFile from file(referenceMap['genomeFile'])
    file genomeIndex from file(referenceMap['genomeIndex'])
    file genomeDict from file(referenceMap['genomeDict'])
    file dbsnp from file(referenceMap['dbsnp'])
    file dbsnpIndex from file(referenceMap['dbsnpIndex'])
    file cosmic from file(referenceMap['cosmic'])
    file cosmicIndex from file(referenceMap['cosmicIndex'])

  output:
    set val("MuTect2"), idPatient, gender, idSampleNormal, idSampleTumor, val("${gen_int}_${idSampleTumor}_vs_${idSampleNormal}"), file("${gen_int}_${idSampleTumor}_vs_${idSampleNormal}.vcf") into mutect2Output

  when: 'MuTect2' in workflowSteps

  // -U ALLOW_SEQ_DICT_INCOMPATIBILITY removed as BAMs generated using the new Picard
  // should be fine

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
  -L \"$genInt\" \
  -XL hs37d5 \
  -XL NC_007605 \
  -o ${gen_int}_${idSampleTumor}_vs_${idSampleNormal}.vcf
  """

}

mutect2Output = mutect2Output.groupTuple(by:[0,1,2,3,4])
verbose ? mutect2Output = mutect2Output.view {"MuTect2 output: $it"} : ''

process RunFreeBayes {
  tag {idPatient + "-" + idSampleTumor + "-" + gen_int}

  input:
    set idPatient, gender, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), genInt, gen_int from bamsFFB
    file genomeFile from file(referenceMap['genomeFile'])

  output:
    set val("FreeBayes"), idPatient, gender, idSampleNormal, idSampleTumor, val("${gen_int}_${idSampleTumor}_vs_${idSampleNormal}"), file("${gen_int}_${idSampleTumor}_vs_${idSampleNormal}.vcf") into freebayesOutput

  when: 'FreeBayes' in workflowSteps

  script:
  """
  freebayes \
    -f $genomeFile \
    --pooled-continuous \
    --pooled-discrete \
    -F 0.03 \
    -C 2 \
    -r \"$genInt\" \
    $bamTumor \
    $bamNormal > ${gen_int}_${idSampleTumor}_vs_${idSampleNormal}.vcf
  """
}

freebayesOutput = freebayesOutput.groupTuple(by:[0,1,2,3,4])
verbose ? freebayesOutput = freebayesOutput.view {"FreeBayes output: $it"} : ''

process RunVardict {
  tag {idPatient + "-" + idSampleTumor + "-" + gen_int}

  input:
    set idPatient, gender, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), genInt, gen_int from bamsFVD
    file genomeFile from file(referenceMap['genomeFile'])
    file genomeIndex from file(referenceMap['genomeIndex'])
    file genomeDict from file(referenceMap['genomeDict'])

  output:
    set val("VarDict"), idPatient, gender, idSampleNormal, idSampleTumor, val("${gen_int}_${idSampleTumor}_vs_${idSampleNormal}"), file("${gen_int}_${idSampleTumor}_vs_${idSampleNormal}.out") into vardictOutput

  when: 'VarDict' in workflowSteps

  script:
  """
  ${referenceMap['vardictHome']}/vardict.pl \
  -G $genomeFile \
  -f 0.01 -N $bamTumor \
  -b "$bamTumor|$bamNormal" \
  -z 1 -F 0x500 \
  -c 1 -S 2 -E 3 -g 4 \
  -R $genInt > ${gen_int}_${idSampleTumor}_vs_${idSampleNormal}.out
  """
}

vardictOutput = vardictOutput.groupTuple(by:[0,1,2,3,4])
verbose ? vardictOutput = vardictOutput.view {"vardictOutput output: $it"} : ''

// we are merging the VCFs that are called separatelly for different intervals
// so we can have a single sorted VCF containing all the calls for a given caller

vcfsToMerge = hcVCF.mix(mutect1Output, mutect2Output, freebayesOutput, vardictOutput)
verbose ? vcfsToMerge = vcfsToMerge.view {"VCFs To be merged: $it"} : ''

process ConcatVCF {
  tag {variantCaller == 'HaplotypeCaller' ? idPatient + "-" + variantCaller + "-" + idSampleNormal : idPatient + "-" + variantCaller + "-" + idSampleNormal + "-" + idSampleTumor}

  publishDir "${directoryMap["$variantCaller"]}", mode: 'copy'

  input:
    set variantCaller, idPatient, gender, idSampleNormal, idSampleTumor, tag, file(vcFiles) from vcfsToMerge
    file genomeFile from file(referenceMap['genomeFile'])
    file genomeDict from file(referenceMap['genomeDict'])
    file genomeIndex from file(referenceMap['genomeIndex'])

  output:
    set variantCaller, idPatient, gender, idSampleNormal, idSampleTumor, file("*.vcf") into vcfConcatenated

  when: 'HaplotypeCaller' in workflowSteps || 'MuTect1' in workflowSteps || 'MuTect2' in workflowSteps || 'FreeBayes' in workflowSteps || 'VarDict' in workflowSteps

  script:
  outputFile = variantCaller == 'HaplotypeCaller' ? "${variantCaller}_${idSampleNormal}.vcf" : "${variantCaller}_${idSampleTumor}_vs_${idSampleNormal}.vcf"
  vcfFiles = vcFiles.collect{" $it"}.join(' ')

  if (variantCaller == 'VarDict')
    """
    for i in $vcFiles ;do
      cat \$i | ${referenceMap['vardictHome']}/VarDict/testsomatic.R >> testsomatic.out
    done
    ${referenceMap['vardictHome']}/VarDict/var2vcf_somatic.pl \
    -f 0.01 \
    -N "${idSampleTumor}_vs_${idSampleNormal}" testsomatic.out > $outputFile
    """

  else if (variantCaller == 'MuTect2' || variantCaller == 'MuTect1' || variantCaller == 'HaplotypeCaller' || variantCaller == 'FreeBayes')
	"""
	# first make a header from one of the VCF intervals
	# get rid of interval information only from the GATK command-line, but leave the rest
	awk '/^#/{print}' `ls *vcf| head -1` | \
	awk '!/GATKCommandLine/{print}/GATKCommandLine/{for(i=1;i<=NF;i++){if(\$i!~/intervals=/ && \$i !~ /out=/){printf("%s ",\$i)}}printf("\\n")}' \
	> header

	## concatenate calls
	rm -rf raw_calls
	for f in *vcf; do
		awk '!/^#/{print}' \$f >> raw_calls
	done
	cat header raw_calls > unsorted.vcf
	java -jar \${PICARD_HOME}/picard.jar SortVcf I=unsorted.vcf O=$outputFile
    rm unsorted.vcf
	"""
}

verbose ? vcfConcatenated = vcfConcatenated.view {"VCF concatenated: $it"} : ''

process RunStrelka {
  tag {idPatient + "-" + idSampleTumor}

  publishDir directoryMap['Strelka'], mode: 'copy'

  input:
    set idPatient, gender, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from bamsForStrelka
    file genomeFile from file(referenceMap['genomeFile'])
    file genomeIndex from file(referenceMap['genomeIndex'])
    file genomeDict from file(referenceMap['genomeDict'])

  output:
    set val("Strelka"), idPatient, gender, idSampleNormal, idSampleTumor, file("*.vcf") into strelkaOutput

  when: 'Strelka' in workflowSteps

  script:
  """
  tumorPath=`readlink $bamTumor`
  normalPath=`readlink $bamNormal`
  genomeFile=`readlink $genomeFile`
  \$STRELKA_INSTALL_DIR/bin/configureStrelkaWorkflow.pl \
  --tumor \$tumorPath \
  --normal \$normalPath \
  --ref \$genomeFile \
  --config \$STRELKA_INSTALL_DIR/etc/strelka_config_bwa_default.ini \
  --output-dir strelka

  cd strelka

  make -j $task.cpus

  cd ..

  mv strelka/results/all.somatic.indels.vcf Strelka_${idSampleTumor}_vs_${idSampleNormal}_all_somatic_indels.vcf
  mv strelka/results/all.somatic.snvs.vcf Strelka_${idSampleTumor}_vs_${idSampleNormal}_all_somatic_snvs.vcf
  mv strelka/results/passed.somatic.indels.vcf Strelka_${idSampleTumor}_vs_${idSampleNormal}_passed_somatic_indels.vcf
  mv strelka/results/passed.somatic.snvs.vcf Strelka_${idSampleTumor}_vs_${idSampleNormal}_passed_somatic_snvs.vcf
  """
}

verbose ? strelkaOutput = strelkaOutput.view {"Strelka output: $it"} : ''

process RunManta {
  tag {idPatient + "-" + idSampleTumor}

  publishDir directoryMap['Manta'], mode: 'copy'

  input:
    set idPatient, gender, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from bamsForManta
    file genomeFile from file(referenceMap['genomeFile'])
    file genomeIndex from file(referenceMap['genomeIndex'])

  output:
    set val("Manta"), idPatient, gender, idSampleNormal, idSampleTumor, file("${idSampleNormal}_${idSampleTumor}.somaticSV.vcf"),file("${idSampleNormal}_${idSampleTumor}.candidateSV.vcf"),file("${idSampleNormal}_${idSampleTumor}.diploidSV.vcf"),file("${idSampleNormal}_${idSampleTumor}.candidateSmallIndels.vcf") into mantaOutput

  when: 'Manta' in workflowSteps

  script:
  """
  set -eo pipefail
  ln -s $bamNormal Normal.bam
  ln -s $bamTumor Tumor.bam
  ln -s $baiNormal Normal.bam.bai
  ln -s $baiTumor Tumor.bam.bai

  configManta.py --normalBam Normal.bam --tumorBam Tumor.bam --reference $genomeFile --runDir MantaDir
  python MantaDir/runWorkflow.py -m local -j $task.cpus
  gunzip -c MantaDir/results/variants/somaticSV.vcf.gz > Manta_${idSampleTumor}_vs_${idSampleNormal}.somaticSV.vcf
  gunzip -c MantaDir/results/variants/candidateSV.vcf.gz > Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSV.vcf
  gunzip -c MantaDir/results/variants/diploidSV.vcf.gz > Manta_${idSampleTumor}_vs_${idSampleNormal}.diploidSV.vcf
  gunzip -c MantaDir/results/variants/candidateSmallIndels.vcf.gz > Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSmallIndels.vcf
  """
}

verbose ? mantaOutput = mantaOutput.view {"Manta output: $it"} : ''

// Run commands and code from Malin Larsson
// Based on Jesper Eisfeldt's code
process RunAlleleCount {
  tag {idPatient + "-" + idSample}

  input:
    set idPatient, gender, status, idSample, file(bam), file(bai) from bamsForAscat
    file genomeFile from file(referenceMap['genomeFile'])
    file genomeIndex from file(referenceMap['genomeIndex'])
    file genomeDict from file(referenceMap['genomeDict'])
    file acLoci from file(referenceMap['acLoci'])

  output:
    set idPatient, gender, status, idSample, file("${idSample}.alleleCount") into alleleCountOutput

  when: 'Ascat' in workflowSteps

  script:
  """
  alleleCounter -l $acLoci -r $genomeFile -b $bam -o ${idSample}.alleleCount;
  """
}

verbose ? alleleCountOutput = alleleCountOutput.view {"alleleCount output: $it"} : ''

alleleCountNormal = Channel.create()
alleleCountTumor = Channel.create()

alleleCountOutput
  .choice(alleleCountTumor, alleleCountNormal) {it[2] =~ /^0$/ ? 1 : 0}

alleleCountOutput = alleleCountNormal.spread(alleleCountTumor)

alleleCountOutput = alleleCountOutput.map {
  idPatientNormal, genderNormal, statusNormal, idSampleNormal, alleleCountNormal, idPatientTumor, genderTumor, statusTumor, idSampleTumor, alleleCountTumor ->
  [idPatientNormal, genderNormal, idSampleNormal, idSampleTumor, alleleCountNormal, alleleCountTumor]
}

verbose ? alleleCountOutput = alleleCountOutput.view {"alleleCount output: $it"} : ''

// R script from Malin Larssons bitbucket repo:
// https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
process RunConvertAlleleCounts {
  tag {idPatient + "-" + idSampleTumor}

  publishDir directoryMap['Ascat'], mode: 'copy'

  input:
    set idPatient, gender, idSampleNormal, idSampleTumor, file(alleleCountNormal), file(alleleCountTumor) from alleleCountOutput

  output:
    set idPatient, gender, idSampleNormal, idSampleTumor, file("${idSampleNormal}.BAF"), file("${idSampleNormal}.LogR"), file("${idSampleTumor}.BAF"), file("${idSampleTumor}.LogR") into convertAlleleCountsOutput

  when: 'Ascat' in workflowSteps

  script:
  """
  convertAlleleCounts.r $idSampleTumor $alleleCountTumor $idSampleNormal $alleleCountNormal $gender
  """
}

// R scripts from Malin Larssons bitbucket repo:
// https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
process RunAscat {
  tag {idPatient + "-" + idSampleTumor}

  publishDir directoryMap['Ascat'], mode: 'copy'

  input:
    set idPatient, gender, idSampleNormal, idSampleTumor, file(bafNormal), file(logrNormal), file(bafTumor), file(logrTumor) from convertAlleleCountsOutput

  output:
    set val("Ascat"), idPatient, gender, idSampleNormal, idSampleTumor, file("${idSampleTumor}.tumour.png"), file("${idSampleTumor}.germline.png"), file("${idSampleTumor}.LogR.PCFed.txt"), file("${idSampleTumor}.BAF.PCFed.txt"), file("${idSampleTumor}.ASPCF.png"), file("${idSampleTumor}.ASCATprofile.png"), file("${idSampleTumor}.aberrationreliability.png"), file("${idSampleTumor}.rawprofile.png"), file("${idSampleTumor}.sunrise.png") into ascatOutput

  when: 'Ascat' in workflowSteps

  script:
  """
  #!/bin/env Rscript
  ##############################################################################
  # Description:                                                               #
  # R-script for converting output from AlleleCount to BAF and LogR values.    #
  #                                                                            #
  # Input:                                                                     #
  # AlleleCounter output file for tumor and normal samples                     #
  # The first line should contain a header describing the data                 #
  # The following columns and headers should be present:                       #
  # CHR    POS     Count_A Count_C Count_G Count_T Good_depth                  #
  #                                                                            #
  # Output:                                                                    #
  # BAF and LogR tables (tab delimited text files)                             #
  ##############################################################################
  source("$baseDir/scripts/ascat.R")
  .libPaths( c( "$baseDir/scripts", .libPaths() ) )
  if(!require(RColorBrewer)){
      source("http://bioconductor.org/biocLite.R")
      biocLite("RColorBrewer", suppressUpdates=TRUE, lib="$baseDir/scripts")
      library(RColorBrewer)
  }

  options(bitmapType='cairo')
  tumorbaf = "$bafTumor"
  tumorlogr = "$logrTumor"
  normalbaf = "$bafNormal"
  normallogr = "$logrNormal"
  #Load the  data
  ascat.bc <- ascat.loadData(Tumor_LogR_file=tumorlogr, Tumor_BAF_file=tumorbaf, Germline_LogR_file=normallogr, Germline_BAF_file=normalbaf)
  #Plot the raw data
  ascat.plotRawData(ascat.bc)
  #Segment the data
  ascat.bc <- ascat.aspcf(ascat.bc)
  #Plot the segmented data
  ascat.plotSegmentedData(ascat.bc)
  #Run ASCAT to fit every tumor to a model, inferring ploidy, normal cell contamination, and discrete copy numbers
  ascat.output <- ascat.runAscat(ascat.bc)
  #str(ascat.output)
  #plot(sort(ascat.output\$aberrantcellfraction))
  #plot(density(ascat.output\$ploidy))
  """
}

verbose ? ascatOutput = ascatOutput.view {"Ascat output: $it"} : ''

// process RunBcftoolsStats {
//   tag {idPatient + "-" + idSample}
//
//   publishDir directoryMap['SamToolsStats'], mode: 'copy'
//
//   input:
//     set idPatient, gender, status, idSample, file(vcf) from ???
//
//   output:
//     file ("${vcf}.samtools.stats.out") into snpeffReport
//
//     when: 'MultiQC' in workflowSteps
//
//     script:
//     """
//     bcfools stats $vcf > ${vcf}.bcf.tools.stats.out
//     """
// }

// process MergeVCF {
//   tag {idPatient + "-" + idSample}
//
//   input:
//     set idPatient, gender, status, idSample, file(vcf) from ???
//
//   output:
//     set idPatient, gender, status, idSample, file ("???") into vcfMerged
//
//     'MuTect1' in workflowSteps || 'MuTect2' in workflowSteps
//
//     script:
//     """
//
//     """
// }

vcfMerged = Channel.create()
vcfNotMerged = Channel.create()

vcfConcatenated
  .choice(vcfMerged, vcfNotMerged) {it[0] == 'MuTect1' ? 0 : 1}

process RunSnpeff {
  tag {variantCaller + "-" + idSampleTumor + "_vs_" + idSampleNormal}

  publishDir directoryMap['snpEff'], mode: 'copy'

  input:
    set variantCaller, idPatient, gender, idSampleNormal, idSampleTumor, file(vcf) from vcfMerged

  output:
    set file("${vcf.baseName}.ann.vcf"), file("${vcf.baseName}_snpEff_genes.txt"), file("${vcf.baseName}_snpEff_summary.html") into snpeffReport

    'MuTect1' in workflowSteps

    script:
    """
    java -Xmx${task.memory.toGiga()}g \
    -jar \$SNPEFF_HOME/snpEff.jar \
    ${params.snpeffDb} \
    -v -cancer \
    ${vcf} \
    > ${vcf.baseName}.ann.vcf

    mv snpEff_genes.txt ${vcf.baseName}_snpEff_genes.txt
    mv snpEff_summary.html ${vcf.baseName}_snpEff_summary.html
    """
}

verbose ? snpeffReport = snpeffReport.view {"snpEff Reports: $it"} : ''

process GenerateMultiQCconfig {
  tag {idPatient}

  publishDir directoryMap['MultiQC'], mode: 'copy'

  input:

  output:
  file("multiqc_config.yaml") into multiQCconfig

  when: 'MultiQC' in workflowSteps

  script:
  """
  touch multiqc_config.yaml
  echo "custom_logo: $baseDir/doc/images/CAW-logo.png" >> multiqc_config.yaml
  echo "custom_logo_url: http://opensource.scilifelab.se/projects/caw" >> multiqc_config.yaml
  echo "custom_logo_title: 'Cancer Analysis Workflow'" >> multiqc_config.yaml
  echo "report_header_info:" >> multiqc_config.yaml
  echo "- CAW version: $version" >> multiqc_config.yaml
  echo "- Contact E-mail: ${params.contactMail}" >> multiqc_config.yaml
  echo "- Command Line: ${workflow.commandLine}" >> multiqc_config.yaml
  echo "- Directory: ${workflow.launchDir}" >> multiqc_config.yaml
  echo "- TSV file: ${tsvFile}" >> multiqc_config.yaml
  echo "- Steps: "${workflowSteps.join(", ")} >> multiqc_config.yaml
  echo "top_modules:" >> multiqc_config.yaml
  echo "- 'fastqc'" >> multiqc_config.yaml
  echo "- 'picard'" >> multiqc_config.yaml
  echo "- 'samtools'" >> multiqc_config.yaml
  echo "- 'snpeff'" >> multiqc_config.yaml
  """
}

verbose ? multiQCconfig = multiQCconfig.view {"MultiQC config file: $it"} : ''

reportsForMultiQC = Channel.fromPath( 'Reports/{FastQC,MarkDuplicates,SamToolsStats}/*' )
  .mix(fastQCreport,markDuplicatesReport,recalibratedBamReport,snpeffReport,multiQCconfig)
  .flatten()
  .unique()
  .toList()

verbose ? reportsForMultiQC = reportsForMultiQC.view {"Reports for MultiQC: $it"} : ''

if (!'MultiQC' in workflowSteps) {reportsForMultiQC.close()}

process RunMultiQC {
  tag {idPatient}

  publishDir directoryMap['MultiQC'], mode: 'copy'

  input:
    file ('*') from reportsForMultiQC

  output:
    set file("*multiqc_report.html"), file("*multiqc_data") into multiQCReport

    when: 'MultiQC' in workflowSteps

  script:
  """
  multiqc -f -v .
  """
}

verbose ? multiQCReport = multiQCReport.view {"MultiQC Report: $it"} : ''


// bamsFVD = bamsFVD.view {"bamsFVD: $it"}

/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def grabGitRevision() { // Borrowed from https://github.com/NBISweden/wgs-structvar
  if (workflow.commitId) { // it's run directly from github
    return workflow.commitId.substring(0,10)
  }
  // Try to find the revision directly from git
  head_pointer_file = file("$baseDir/.git/HEAD")
  if (!head_pointer_file.exists()) {
    return ''
  }
  ref = head_pointer_file.newReader().readLine().tokenize()[1]
  ref_file = file("$baseDir/.git/$ref")
  if (!ref_file.exists()) {
    return ''
  }
  revision = ref_file.newReader().readLine()
  return revision.substring(0,10)
}

def checkUppmaxProject() {
  if ((workflow.profile == 'standard' || workflow.profile == 'interactive') && !params.project) {
    return false
  } else {
    return true
  }
}

def defineReferenceMap() {
  return [
    'acLoci'      : params.acLoci,      // loci file for ascat
    'dbsnp'       : params.dbsnp,       // dbSNP
    'cosmic'      : params.cosmic,      // cosmic vcf file with VCF4.1 header
    'cosmicIndex' : params.cosmicIndex, // cosmic vcf file index
    'dbsnpIndex'  : params.dbsnpIndex,  // dbSNP index
    'genomeAmb'   : params.genomeAmb,   // BWA indexes
    'genomeAnn'   : params.genomeAnn,   // BWA indexes
    'genomeBwt'   : params.genomeBwt,   // BWA indexes
    'genomeDict'  : params.genomeDict,  // genome reference dictionary
    'genomeFile'  : params.genome,      // genome reference
    'genomeIndex' : params.genomeIndex, // genome reference index
    'genomePac'   : params.genomePac,   // BWA indexes
    'genomeSa'    : params.genomeSa,    // BWA indexes
    'intervals'   : params.intervals,   // intervals file for spread-and-gather processes (usually chromosome chunks at centromeres)
    'kgIndels'    : params.kgIndels,    // 1000 Genomes SNPs
    'kgIndex'     : params.kgIndex,     // 1000 Genomes SNPs index
    'millsIndels' : params.millsIndels, // Mill's Golden set of SNPs
    'millsIndex'  : params.millsIndex,  // Mill's Golden set index
    'vardictHome' : params.vardictHome  // path to VarDict
  ]
}

def defineDirectoryMap() {
  return [
    'nonRealigned'     : 'Preprocessing/NonRealigned',
    'nonRecalibrated'  : 'Preprocessing/NonRecalibrated',
    'recalibrated'     : 'Preprocessing/Recalibrated',
    'FastQC'           : 'Reports/FastQC',
    'MarkDuplicatesQC' : 'Reports/MarkDuplicates',
    'MultiQC'          : 'Reports/MultiQC',
    'SamToolsStats'    : 'Reports/SamToolsStats',
    'Ascat'            : 'VariantCalling/Ascat',
    'FreeBayes'        : 'VariantCalling/FreeBayes',
    'HaplotypeCaller'  : 'VariantCalling/HaplotypeCaller',
    'Manta'            : 'VariantCalling/Manta',
    'MuTect1'          : 'VariantCalling/MuTect1',
    'MuTect2'          : 'VariantCalling/MuTect2',
    'Strelka'          : 'VariantCalling/Strelka',
    'VarDict'          : 'VariantCalling/VarDict',
    'snpEff'           : 'Annotation/snpEff'
  ]
}

def defineStepList() {
  return [
    'Ascat',
    'FreeBayes',
    'HaplotypeCaller',
    'Manta',
    'MultiQC',
    'MuTect1',
    'MuTect2',
    'preprocessing',
    'realign',
    'recalibrate',
    'skipPreprocessing',
    'Strelka',
    'VarDict'
  ]
}

def checkReferenceMap(referenceMap) {
  // Loop through all the references files to check their existence
  referenceDefined = true
  referenceMap.each{
    referenceFile, fileToCheck ->
    test = checkRefExistence(referenceFile, fileToCheck)
    !(test) ? referenceDefined = false : ''
  }
  return referenceDefined ? true : false
}

def checkStepList(stepsList, realStepsList) {
  // Loop through all the possible steps check their existence and spelling
  stepCorrect = true
  stepsList.each{
    test = checkStepExistence(it, realStepsList)
    !(test) ? stepCorrect = false : ''
  }
  return stepCorrect ? true : false
}

def checkRefExistence(referenceFile, fileToCheck) {
  // Check file existence
  try {assert file(fileToCheck).exists()}
  catch (AssertionError ae) {
    log.info  "Missing references: $referenceFile $fileToCheck"
    return false
  }
  return true
}

def checkStepExistence(step, list) {
  // Check step existence
  try {assert list.contains(step)}
  catch (AssertionError ae) {
    println("Unknown parameter: $step")
    return false
  }
  return true
}

def checkFile(it) {
  // Check file existence
  try {assert file(it).exists()}
  catch (AssertionError ae) {
    exit 1, "Missing file in TSV file: $it, see --help for more information"
  }
  return it
}

def checkTSV(it,number) {
  // Check if TSV has the correct number of item in row
  try {assert it.size() == number}
  catch (AssertionError ae) {
    exit 1, "Malformed row in TSV file: $it, see --help for more information"
  }
  return it
}

def checkStatus(it) {
  // Check if Status is correct
  try {assert it == "0" || it == "1"}
  catch (AssertionError ae) {
    exit 1, "Status is not recognized in TSV file: $it, see --help for more information"
  }
  return it
}

def extractFastqFiles(tsvFile) {
  // Channeling the TSV file containing FASTQ.
  // Format is: "subject gender status sample lane fastq1 fastq2"
  fastqFiles = Channel
    .from(tsvFile.readLines())
    .map{line ->
      list       = checkTSV(line.split(),7)
      idPatient  = list[0]
      gender     = list[1]
      status     = checkStatus(list[2])
      idSample   = list[3]
      idRun      = list[4]

      // When testing workflow from github, paths to FASTQ files start from workflow.projectDir and not workflow.launchDir
      fastqFile1 = workflow.commitId && params.test ? checkFile("$workflow.projectDir/${list[5]}") : checkFile("${list[5]}")
      fastqFile2 = workflow.commitId && params.test ? checkFile("$workflow.projectDir/${list[6]}") : checkFile("${list[6]}")

      [idPatient, gender, status, idSample, idRun, fastqFile1, fastqFile2]
    }
  return fastqFiles
}

def extractBamFiles(tsvFile) {
  // Channeling the TSV file containing BAM.
  // Format is: "subject gender status sample bam bai"
  bamFiles = Channel
    .from(tsvFile.readLines())
    .map{line ->
      list      = checkTSV(line.split(),6)
      idPatient = list[0]
      gender    = list[1]
      status    = checkStatus(list[2])
      idSample  = list[3]
      bamFile   = checkFile(list[4])
      baiFile   = checkFile(list[5])

      [ idPatient, gender, status, idSample, bamFile, baiFile ]
    }
  return bamFiles
}

def extractRecalibrationTables(tsvFile) {
  // Channeling the TSV file containing Recalibration Tables.
  // Format is: "subject gender status sample bam bai recalTables"
  bamFiles = Channel
    .from(tsvFile.readLines())
    .map{line ->
      list       = checkTSV(line.split(),7)
      idPatient  = list[0]
      gender     = list[1]
      status     = checkStatus(list[2])
      idSample   = list[3]
      bamFile    = checkFile(list[4])
      baiFile    = checkFile(list[5])
      recalTable = checkFile(list[6])

      [ idPatient, gender, status, idSample, bamFile, baiFile, recalTable ]
    }
  return bamFiles
}

def checkSteps(workflowSteps) {
  result = 0

  if ('preprocessing' in workflowSteps) {result++}
  if ('recalibrate' in workflowSteps) {result++}
  if ('realign' in workflowSteps) {result++}
  if ('skipPreprocessing' in workflowSteps) {result++}
  if (result == 1 ) {
    return true
  } else {
    return false
  }
}

def retreiveStatus(bamChannel) {
  bamChannel = bamChannel.map {
    idPatient, gender, bam, bai ->
    tag = bam.baseName.tokenize('.')[0]
    status   = tag[-1..-1]
    idSample = tag.take(tag.length()-2)
    [idPatient, gender, status, idSample, bam, bai]
  }
  return bamChannel
}

def generateIntervalsForVC(bams, gI) {
  (bams, bamsForVC) = bams.into(2)
  (gI, vcIntervals) = gI.into(2)
  bamsForVC = bamsForVC.spread(vcIntervals)
  return [bamsForVC, bams, gI]
}

def help_message(version, revision) { // Display help message
  log.info "CANCER ANALYSIS WORKFLOW ~ $version - revision: $revision"
  log.info "    Usage:"
  log.info "       nextflow run SciLifeLab/CAW --sample <sample.tsv> [--steps STEP[,STEP]]"
  log.info "    --steps"
  log.info "       Option to configure which processes to use in the workflow."
  log.info "         Different steps to be separated by commas."
  log.info "       Possible values are:"
  log.info "         preprocessing (default, will start workflow with FASTQ files)"
  log.info "         recalibrate (will start workflow with non-recalibrated BAM files)"
  log.info "         skipPreprocessing (will start workflow with recalibrated BAM files)"
  log.info "         MuTect1 (use MuTect1 for VC)"
  log.info "         MuTect2 (use MuTect2 for VC)"
  log.info "         FreeBayes (use FreeBayes for VC)"
  log.info "         VarDict (use VarDict for VC)"
  log.info "         Strelka (use Strelka for VC)"
  log.info "         HaplotypeCaller (use HaplotypeCaller for normal bams VC)"
  log.info "         Manta (use Manta for SV)"
  log.info "         Ascat (use Ascat for CNV)"
  log.info "    --help"
  log.info "       you're reading it"
  log.info "    --verbose"
  log.info "       Adds more verbosity to workflow"
  log.info "    --version"
  log.info "       displays version number"
  log.info "    Test:"
  log.info "      to test CAW on smaller dataset, enter one of the following command"
  log.info "    nextflow run SciLifeLab/CAW --test"
  log.info "       Test `preprocessing` on test tiny set"
  log.info "    nextflow run SciLifeLab/CAW --testRealign"
  log.info "       Test `realign` on test tiny set"
  log.info "       Need to do `nextflow run SciLifeLab/CAW --test` before"
  log.info "    nextflow run SciLifeLab/CAW --testCoreVC"
  log.info "       Test `preprocessing`, `MuTect1`, `Strelka` and `HaplotypeCaller` on test tiny set"
  log.info "       Need to do `nextflow run SciLifeLab/CAW --test` before"
  log.info "    nextflow run SciLifeLab/CAW --testSideVC"
  log.info "       Test `skipPreprocessing`, `Ascat`, `Manta` and `HaplotypeCaller` on test downSampled set"
}

def start_message(version, revision) { // Display start message
  log.info "CANCER ANALYSIS WORKFLOW ~ $version - revision: $revision"
  log.info "Command Line: $workflow.commandLine"
  log.info "Project Dir : $workflow.projectDir"
  log.info "Launch Dir  : $workflow.launchDir"
  log.info "Work Dir    : $workflow.workDir"
  log.info "Steps       : " + workflowSteps.join(', ')
}

def version_message(version, revision) { // Display version message
  log.info "CANCER ANALYSIS WORKFLOW"
  log.info "  version   : $version"
  log.info workflow.commitId ? "Git info    : $workflow.repository - $workflow.revision [$workflow.commitId]" : "  revision  : $revision"
}

workflow.onComplete { // Display complete message
  log.info "N E X T F L O W ~ $workflow.nextflow.version - $workflow.nextflow.build"
  log.info "CANCER ANALYSIS WORKFLOW ~ $version - revision: $revision"
  log.info "Command Line: $workflow.commandLine"
  log.info "Project Dir : $workflow.projectDir"
  log.info "Launch Dir  : $workflow.launchDir"
  log.info "Work Dir    : $workflow.workDir"
  log.info "TSV file    : $tsvFile"
  log.info "Steps       : " + workflowSteps.join(", ")
  log.info "Completed at: $workflow.complete"
  log.info "Duration    : $workflow.duration"
  log.info "Success     : $workflow.success"
  log.info "Exit status : $workflow.exitStatus"
  log.info "Error report: " + (workflow.errorReport ?: '-')
}

workflow.onError { // Display error message
  log.info "N E X T F L O W ~ version $workflow.nextflow.version [$workflow.nextflow.build]"
  log.info workflow.commitId ? "CANCER ANALYSIS WORKFLOW ~ $version - $workflow.revision [$workflow.commitId]" : "CANCER ANALYSIS WORKFLOW ~ $version - revision: $revision"
  log.info "Workflow execution stopped with the following message: " + workflow.errorMessage
}
