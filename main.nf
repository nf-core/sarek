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
 https://github.com/SciLifeLab/CAW
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
 - RecalibrateBam - Recalibreate Bam
 - RunHaplotypecaller - Run HaplotypeCaller for GermLine Variant Calling (Parrallelized processes)
 - RunMutect1 - Run MuTect1 for Variant Calling (Parrallelized processes)
 - RunMutect2 - Run MuTect2 for Variant Calling (Parrallelized processes)
 - RunFreeBayes - Run FreeBayes for Variant Calling (Parrallelized processes)
 - RunVardict - Run VarDict for Variant Calling (Parrallelized processes)
 - ConcatVCF - Merge results from HaplotypeCaller, MuTect1, MuTect2 and VarDict
 - RunStrelka - Run Strelka for Variant Calling
 - RunManta - Run Manta for Structural Variant Calling
 - RunAlleleCount - Run AlleleCount to prepare for Ascat
 - RunConvertAlleleCounts - Run convertAlleleCounts to prepare for Ascat
 - RunAscat - Run Ascat for CNV
 - RunMultiQC - Run MultiQC for report and QC
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

revision = grabGitRevision() ?: ''
version = 'v0.9.9'
verbose = false
testFile = ''
testSteps = []
workflowSteps = []

if (!checkUppmaxProject()) {exit 1, 'No UPPMAX project ID found! Use --project <UPPMAX Project ID>'}

switch (params) {
  case {params.help} :
    help_message(version, revision)
    exit 1

  case {params.version} :
    version_message(version, revision)
    exit 1

  // Arguments handling: Getting list of steps from comma-separated strings to choose which processes to run or to skip
  // Borrowed from https://github.com/guigolab/grape-nf and https://github.com/NBISweden/wgs-structvar
  case {params.steps} :
    workflowSteps = params.steps.split(',').collect {it.trim()}
    break
}

if (params.verbose) {verbose = true}

referenceMap = defineReferenceMap()
directoryMap = defineDirectoryMap()
stepList = defineStepList()

if (!checkReferenceMap(referenceMap)) {exit 1, 'Missing Reference file(s), see --help for more information'}
if (!checkStepList(workflowSteps,stepList)) {exit 1, 'Unknown step(s), see --help for more information'}

if (params.test) {
  test = true
  testFile = file("$workflow.projectDir/data/tsv/tiny.tsv")
  workflowSteps = ['preprocessing']
  referenceMap.put("intervals", "$workflow.projectDir/repeats/tiny.list")
} else if (params.testRealign) {
  test = true
  testFile = file("$workflow.launchDir/${directoryMap['nonRealigned']}/nonRealigned.tsv")
  workflowSteps = ['realign']
  referenceMap.put("intervals", "$workflow.projectDir/repeats/tiny.list")
} else if (params.testCoreVC) {
  test = true
  testFile = file("$workflow.launchDir/${directoryMap['recalibrated']}/recalibrated.tsv")
  workflowSteps = ['skipPreprocessing', 'MuTect1', 'Strelka', 'HaplotypeCaller']
  referenceMap.put("intervals", "$workflow.projectDir/repeats/tiny.list")
} else if (params.testSideVC) {
  test = true
  testFile = file("$workflow.projectDir/data/tsv/G15511-recalibrated.tsv")
  workflowSteps = ['skipPreprocessing', 'Ascat', 'Manta', 'HaplotypeCaller']
} else {test = false}

if (('preprocessing' in workflowSteps && ('realign' in workflowSteps || 'skipPreprocessing' in workflowSteps)) || ('realign' in workflowSteps && 'skipPreprocessing' in workflowSteps)) {
  exit 1, 'Please choose only one step between preprocessing, realign and skipPreprocessing, see --help for more information'
}
if (!('preprocessing' in workflowSteps || 'realign' in workflowSteps || 'skipPreprocessing' in workflowSteps)) {
  exit 1, 'Please choose one step between preprocessing, realign and skipPreprocessing, see --help for more information'
}

/*
 * Extract and verify content of TSV file
 */

if ((!params.sample) && !(test)) {exit 1, 'Missing TSV file, see --help for more information'}
tsvFile = (!(test) ? file(params.sample) : testFile)

fastqFiles = Channel.create()
fastqFilesforFastQC = Channel.create()

if ('preprocessing' in workflowSteps) {
  fastqFiles = extractFastqFiles(tsvFile)
  if (verbose) {fastqFiles = fastqFiles.view {"FASTQ files and IDs to process: $it"}}
} else if ('realign' in workflowSteps || 'skipPreprocessing' in workflowSteps) {
  bamFiles = extractBamFiles(tsvFile)
  if (verbose) {bamFiles = bamFiles.view {"Bam files and IDs to process: $it"}}
  fastqFiles.close()
}

start_message(version, revision)

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

if ('preprocessing' in workflowSteps && 'MultiQC' in workflowSteps) {
  (fastqFiles, fastqFilesforFastQC) = fastqFiles.into(2)
}

process RunFastQC {
  tag {idRun}

  input:
    set idPatient, gender, status, idSample, idRun, file(fastqFile1), file(fastqFile2) from fastqFilesforFastQC

  output:
    file "*_fastqc.{zip,html}" into fastQCreport

  when: 'preprocessing' in workflowSteps && 'MultiQC' in workflowStepsworkflowSteps

  script:
  """
  fastqc -q $fastqFile1 $fastqFile2
  """
}

if ('preprocessing' in workflowSteps && 'MultiQC' in workflowSteps) {
  if (verbose) {fastQCreport = fastQCreport.view {"FastQC report: $it"}}
}

process MapReads {
  tag {idRun}

  input:
    set idPatient, gender, status, idSample, idRun, file(fastqFile1), file(fastqFile2) from fastqFiles

  output:
    set idPatient, gender, status, idSample, idRun, file("${idRun}.bam") into bam

  when: 'preprocessing' in workflowSteps

  script:
  readGroup="@RG\\tID:$idRun\\tSM:$idSample\\tLB:$idSample\\tPL:illumina"
  """
  set -eo pipefail
  bwa mem -R \"$readGroup\" -B 3 -t $task.cpus -M ${referenceMap['genomeFile']} $fastqFile1 $fastqFile2 | \
  samtools view -bS -t ${referenceMap['genomeIndex']} - | \
  samtools sort - > ${idRun}.bam
  """
}

singleBam = Channel.create()
groupedBam = Channel.create()

if ('preprocessing' in workflowSteps) {
  if (verbose) {bam = bam.view {"BAM file before sorting into group or single: $it"}}

  /*
   * Sort bam whether they are standalone or should be merged
   * Borrowed code from https://github.com/guigolab/chip-nf
   */

  bam.groupTuple(by:[0,1,2,3])
    .choice(singleBam, groupedBam) {it[4].size() > 1 ? 1 : 0}
  singleBam = singleBam.map {
    idPatient, gender, status, idSample, idRun, bam ->
    [idPatient, gender, status, idSample, bam]
  }
  if (verbose) {groupedBam = groupedBam.view {"Grouped BAMs before merge: $it"}}
} else {
  singleBam.close()
  groupedBam.close()
}

process MergeBams {
  tag {idSample}

  input:
    set idPatient, gender, status, idSample, idRun, file(bam) from groupedBam

  output:
    set idPatient, gender, status, idSample, file("${idSample}.bam") into mergedBam

  when: 'preprocessing' in workflowSteps

  script:
  """
  samtools merge ${idSample}.bam $bam
  """
}

if ('preprocessing' in workflowSteps) {
  /*
   * Gather all bams into a single channel
   */

  if (verbose) {
    singleBam = singleBam.view {"Single BAM: $it"}
    mergedBam = mergedBam.view {"Merged BAM: $it"}
  }
  mergedBam = mergedBam.mix(singleBam)
  if (verbose) {mergedBam = mergedBam.view {"BAM for MarkDuplicates: $it"}}
}
process MarkDuplicates {
  tag {idSample}

  publishDir directoryMap['nonRealigned'], mode: 'copy'

  input:
    set idPatient, gender, status, idSample, file(bam) from mergedBam

  output:
    set idPatient, gender, val("${idSample}_${status}"), file("${idSample}_${status}.md.bam"), file("${idSample}_${status}.md.bai") into duplicates
    set idPatient, gender, status, idSample, val("${idSample}_${status}.md.bam"), val("${idSample}_${status}.md.bai") into markDuplicatesTSV

  when: 'preprocessing' in workflowSteps

  script:
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar ${referenceMap['picardHome']}/MarkDuplicates.jar \
  INPUT=${bam} \
  METRICS_FILE=${bam}.metrics \
  TMP_DIR=. \
  ASSUME_SORTED=true \
  VALIDATION_STRINGENCY=LENIENT \
  CREATE_INDEX=TRUE \
  OUTPUT=${idSample}_${status}.md.bam
  """
}

markDuplicatesTSV.map { idPatient, gender, status, idSample, bam, bai ->
  "$idPatient\t$gender\t$status\t$idSample\t${directoryMap['nonRealigned']}/$bam\t${directoryMap['nonRealigned']}/$bai\n"
}.collectFile( name: 'nonRealigned.tsv', sort: true, storeDir: directoryMap['nonRealigned'])

duplicatesInterval = Channel.create()
duplicatesRealign  = Channel.create()

if ('preprocessing' in workflowSteps) {
  duplicatesGrouped  = Channel.create()
  // create intervals for realignement using both tumor+normal as input
  // group the marked duplicates Bams for intervals and realign by overall subject/patient id (idPatient)
  duplicatesGrouped = duplicates.groupTuple(by:[0,1])
} else if ('realign' in workflowSteps) {
  duplicatesGrouped = Channel.create()
  duplicatesGrouped = bamFiles.map{
    idPatient, gender, status, idSample, bam, bai ->
    [idPatient, gender, "${idSample}_${status}", bam, bai]
  }.groupTuple(by:[0,1])
}

if ('preprocessing' in workflowSteps || 'realign' in workflowSteps) {
  // The duplicatesGrouped channel is duplicated, one copy goes to the CreateIntervals process
  // and the other to the IndelRealigner process
  duplicatesGrouped.into(duplicatesInterval, duplicatesRealign)
  if (verbose) {
    duplicatesInterval = duplicatesInterval.view {"BAMs for RealignerTargetCreator grouped by patient ID: $it"}
    duplicatesRealign  = duplicatesRealign.view  {"BAMs for IndelRealigner grouped by patient ID: $it"}
  }
} else {
  duplicatesInterval.close()
  duplicatesRealign.close()
}

// VCF indexes are added so they will be linked, and not re-created on the fly
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
    file("${idPatient}.intervals") into intervals

  when: 'preprocessing' in workflowSteps || 'realign' in workflowSteps

  script:
  bams = bam.collect{"-I $it"}.join(' ')
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar ${referenceMap['gatkHome']}/GenomeAnalysisTK.jar \
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

if ('preprocessing' in workflowSteps || 'realign' in workflowSteps) {
  if (verbose) {intervals = intervals.view {"Intervals passed to Realign: $it"}}
}

// use nWayOut to split into T/N pair again
// VCF indexes are added so they will be linked, and not re-created on the fly
process RealignBams {
  tag {idPatient}

  input:
    set idPatient, gender, idSample_status, file(bam), file(bai) from duplicatesRealign
    file genomeFile from file(referenceMap['genomeFile'])
    file genomeIndex from file(referenceMap['genomeIndex'])
    file genomeDict from file(referenceMap['genomeDict'])
    file kgIndels from file(referenceMap['kgIndels'])
    file kgIndex from file(referenceMap['kgIndex'])
    file millsIndels from file(referenceMap['millsIndels'])
    file millsIndex from file(referenceMap['millsIndex'])
    file intervals from intervals

  output:
    val(idPatient) into tempIdPatient
    val(gender) into tempGender
    val(idSample_status) into tempSamples_status
    file("*.real.bam") into tempBams
    file("*.real.bai") into tempBais

  when: 'preprocessing' in workflowSteps || 'realign' in workflowSteps

  script:
  bams = bam.collect{"-I $it"}.join(' ')
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar ${referenceMap['gatkHome']}/GenomeAnalysisTK.jar \
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

realignedBam = Channel.create()

if ('preprocessing' in workflowSteps || 'realign' in workflowSteps) {
  // If I make a set out of this process I got a list of lists, which cannot be iterate via a single process
  // So I need to transform this output into a channel that can be iterated on.
  // I also had problems with the set that wasn't synchronised, and I got wrongly associated files
  // So what I decided was to separate all the different output
  // We're getting from the Realign process 5 channels (patient, gender, samples bams and bais)
  // So I flatten, sort, and reflatten the samples and the files (bam and bai) channels
  // to get them in the same order (the name of the bam and bai files are based on the sample, so if we sort them they all have the same order ;-))
  // And put them back together, and add the ID patient and the gender in the realignedBam channel
  realignedBam = tempIdPatient.spread(
    tempGender.spread(
      tempSamples_status.flatten().toSortedList().flatten().merge(
        tempBams.flatten().toSortedList().flatten(),
        tempBais.flatten().toSortedList().flatten()
      ) {sample, bam, bai -> [sample, bam, bai]}
    )
  )
  //Retrieving the status from the idSample
  realignedBam = retreiveStatus(realignedBam)
  if (verbose) {realignedBam = realignedBam.view {"realignedBam to BaseRecalibrator: $it"}}
} else {
  realignedBam.close()
}

process CreateRecalibrationTable {
  tag {idSample}

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

  when: 'preprocessing' in workflowSteps || 'realign' in workflowSteps

  script:
  """
  java -Xmx${task.memory.toGiga()}g \
  -Djava.io.tmpdir="/tmp" \
  -jar ${referenceMap['gatkHome']}/GenomeAnalysisTK.jar \
  -T BaseRecalibrator \
  -R $genomeFile \
  -I $bam \
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

if ('preprocessing' in workflowSteps || 'realign' in workflowSteps) {
  if (verbose) {recalibrationTable = recalibrationTable.view {"Base recalibrated table for recalibration: $it"}}
}

process RecalibrateBam {
  tag {idSample}

  publishDir directoryMap['recalibrated'], mode: 'copy'

  input:
    set idPatient, gender, status, idSample, file(bam), file(bai), recalibrationReport from recalibrationTable
    file genomeFile from file(referenceMap['genomeFile'])
    file genomeIndex from file(referenceMap['genomeIndex'])
    file genomeDict from file(referenceMap['genomeDict'])

  output:
    set idPatient, gender, status, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bai") into recalibratedBam
    set idPatient, gender, status, idSample, val("${idSample}.recal.bam"), val("${idSample}.recal.bai") into recalibratedBamTSV

  when: 'preprocessing' in workflowSteps || 'realign' in workflowSteps

  // TODO: ditto as at the previous BaseRecalibrator step, consider using -nct 4
  script:
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar ${referenceMap['gatkHome']}/GenomeAnalysisTK.jar \
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

recalibratedBamTSV.map { idPatient, gender, status, idSample, bam, bai ->
  "$idPatient\t$gender\t$status\t$idSample\t${directoryMap['recalibrated']}/$bam\t${directoryMap['recalibrated']}/$bai\n"
}.collectFile( name: 'recalibrated.tsv', sort: true, storeDir: directoryMap['recalibrated'])

if ('skipPreprocessing' in workflowSteps) {
  recalibratedBam = bamFiles
}

if (verbose) {recalibratedBam = recalibratedBam.view {"Recalibrated Bam for variant Calling: $it"}}

// Here we have a recalibrated bam set, but we need to separate the bam files based on patient status.
// The sample tsv config file which is formatted like: "subject status sample lane fastq1 fastq2"
// cf fastqFiles channel, I decided just to add _status to the sample name to have less changes to do.
// And so I'm sorting the channel if the sample match _0, then it's a normal sample, otherwise tumor.
// Then spread normal over tumor to get each possibilities
// ie. normal vs tumor1, normal vs tumor2, normal vs tumor3
// then copy this channel into channels for each variant calling
// I guess it will still work even if we have multiple normal samples

bamsTumor = Channel.create()
bamsNormal = Channel.create()
bamsAll = Channel.create()
vcfsToMerge = Channel.create()

bamsFHC = Channel.create()  // HaplotypeCaller
bamsFMT1 = Channel.create() // MuTect1
bamsFMT2 = Channel.create() // MuTect2
bamsFFB = Channel.create()  // FreeBayes
bamsFVD = Channel.create() // VarDict

hcVCF = Channel.create()        // HaplotypeCaller
mutect1VCF = Channel.create()   // MuTect1
mutect2VCF = Channel.create()   // MuTect2
freebayesVCF = Channel.create() // FreeBayes
vardictVCF = Channel.create()   // VarDict

bamsForStrelka = Channel.create()
bamsForManta = Channel.create()
bamsForAscat = Channel.create()

// separate recalibrateBams by status

recalibratedBam
  .choice(bamsTumor, bamsNormal) {it[2] =~ /^0$/ ? 1 : 0}

// Removing status because not relevant anymore
bamsTumor = bamsTumor.map { idPatient, gender, status, idSample, bam, bai -> [idPatient, gender, idSample, bam, bai] }
bamsNormal = bamsNormal.map { idPatient, gender, status, idSample, bam, bai -> [idPatient, gender, idSample, bam, bai] }

if (verbose) {
  bamsTumor = bamsTumor.view {"Tumor Bam for variant Calling: $it"}
  bamsNormal = bamsNormal.view {"Normal Bam for variant Calling: $it"}
}

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

if ('HaplotypeCaller' in workflowSteps) {
  bamsFHCTemp = Channel.create()
  (bamsFHC, bamsNormal, gI) = generateIntervalsForVC(bamsNormal, gI)
  (bamsFHCTemp, bamsTumor, gI) = generateIntervalsForVC(bamsTumor, gI)
  bamsFHC = bamsFHC.mix(bamsFHCTemp)
  if (verbose) {bamsFHC = bamsFHC.view {"Bams with Intervals for HaplotypeCaller: $it"}}
} else {
  bamsFHC.close()
  hcVCF.close()
}

bamsAll = bamsNormal.spread(bamsTumor)

bamsAll = bamsAll.map { // Since idPatientNormal and idPatientTumor are the same, it's removed from bamsAll Channel (same for genderNormal)
  // /!\ It is assumed that every sample are from the same patient
  idPatientNormal, genderNormal, idSampleNormal, bamNormal, baiNormal, idPatientTumor, genderTumor, idSampleTumor, bamTumor, baiTumor ->
  [idPatientNormal, genderNormal, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor]
}

if (verbose) {bamsAll = bamsAll.view {"Mapped Recalibrated Bam for variant Calling: $it"}}

if ('MuTect1' in workflowSteps) {
  (bamsFMT1, bamsAll, gI) = generateIntervalsForVC(bamsAll, gI)
  if (verbose) {bamsFMT1 = bamsFMT1.view {"Bams with Intervals for MuTect1: $it"}}
} else {
  bamsFMT1.close()
  mutect1VCF.close()
}

if ('MuTect2' in workflowSteps) {
  (bamsFMT2, bamsAll, gI) = generateIntervalsForVC(bamsAll, gI)
  if (verbose) {bamsFMT2 = bamsFMT2.view {"Bams with Intervals for MuTect2: $it"}}
} else {
  bamsFMT2.close()
  mutect2VCF.close()
}

if ('FreeBayes' in workflowSteps) {
  (bamsFFB, bamsAll, gI) = generateIntervalsForVC(bamsAll, gI)
  if (verbose) {bamsFFB = bamsFFB.view {"Bams with Intervals for FreeBayes: $it"}}
} else {
  bamsFFB.close()
  freebayesVCF.close()
}

if ('VarDict' in workflowSteps) {
  (bamsFVD, bamsAll, gI) = generateIntervalsForVC(bamsAll, gI)
  if (verbose) {bamsFVD = bamsFVD.view {"Bams with Intervals for VarDict: $it"}}
} else {
  bamsFVD.close()
  vardictVCF.close()
}

if ('Strelka' in workflowSteps) {
  (bamsAll, bamsForStrelka) = bamsAll.into(2)
  if (verbose) {bamsForStrelka = bamsForStrelka.view {"Bams for Strelka: $it"}}
} else {
  bamsForStrelka.close()
}

if ('Manta' in workflowSteps) {
  (bamsAll, bamsForManta) = bamsAll.into(2)
  if (verbose) {bamsForManta = bamsForManta.view {"Bams for Manta: $it"}}
} else {
  bamsForManta.close()
}

if ('Ascat' in workflowSteps) {
  (bamsAll, bamsForAscat) = bamsAll.into(2)
  if (verbose) {bamsForAscat = bamsForAscat.view {"Bams for Ascat: $it"}}
} else {
  bamsForAscat.close()
}

process RunHaplotypecaller {
  tag {idSample + "-" + gen_int}

  input:
    set idPatient, gender, idSample, file(bam), file(bai), genInt, gen_int from bamsFHC //Are these values `ped to bamNormal already?
    file genomeFile from file(referenceMap['genomeFile'])
    file genomeIndex from file(referenceMap['genomeIndex'])
    file genomeDict from file(referenceMap['genomeDict'])
    file dbsnp from file(referenceMap['dbsnp'])
    file dbsnpIndex from file(referenceMap['dbsnpIndex'])

  output:
    set val("HaplotypeCaller"), idPatient, gender, idSample, val("${gen_int}_${idSample}"), file("${gen_int}_${idSample}.vcf") into haplotypecallerOutput

  when: 'HaplotypeCaller' in workflowSteps

  //parellelization information: "Many users have reported issues running HaplotypeCaller with the -nct argument, so we recommend using Queue to parallelize HaplotypeCaller instead of multithreading." However, it can take the -nct argument.
  script:
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar ${referenceMap['gatkHome']}/GenomeAnalysisTK.jar \
  -T HaplotypeCaller \
  -R $genomeFile \
  --dbsnp $dbsnp \
  -I $bam \
  -L \"$genInt\" \
  -XL hs37d5 \
  -XL NC_007605 \
  -o ${gen_int}_${idSample}.vcf
  """
}

if ('HaplotypeCaller' in workflowSteps) {
  if (verbose) {haplotypecallerOutput = haplotypecallerOutput.view {"HaplotypeCaller output: $it"}}
  hcVCF = haplotypecallerOutput.map {
    variantCaller, idPatient, gender, idSample, tag, vcfFile ->
    [variantCaller, idPatient, gender, idSample, idSample, tag, vcfFile]
  }.groupTuple(by:[0,1,2,3,4])
  if (verbose) {hcVCF = hcVCF.view {"hcVCF: $it" } }
}

process RunMutect1 {
  tag {idSampleTumor + "-" + gen_int}

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
    set val("MuTect1"), idPatient, gender, idSampleNormal, idSampleTumor, val("${gen_int}_${idSampleNormal}_${idSampleTumor}"), file("${gen_int}_${idSampleNormal}_${idSampleTumor}.vcf") into mutect1Output

  when: 'MuTect1' in workflowSteps

  script:
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar ${referenceMap['mutect1Home']}/muTect.jar \
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
  --out ${gen_int}_${idSampleNormal}_${idSampleTumor}.call_stats.out \
  --vcf ${gen_int}_${idSampleNormal}_${idSampleTumor}.vcf
  """
}

if ('MuTect1' in workflowSteps) {
  if (verbose) {mutect1Output = mutect1Output.view {"MuTect1 output: $it"}}
  mutect1VCF = mutect1Output.groupTuple(by:[0,1,2,3,4])
  if (verbose) {mutect1VCF = mutect1VCF.view {"mutect1VCF: $it"}}
}

process RunMutect2 {
  tag {idSampleTumor + "-" + gen_int}

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
    set val("MuTect2"), idPatient, gender, idSampleNormal, idSampleTumor, val("${gen_int}_${idSampleNormal}_${idSampleTumor}"), file("${gen_int}_${idSampleNormal}_${idSampleTumor}.vcf") into mutect2Output

  // we are using MuTect2 shipped in GATK v3.6
  // TODO: the  "-U ALLOW_SEQ_DICT_INCOMPATIBILITY " flag is actually masking a bug in older Picard versions. Using the latest Picard tool
  // this bug should go away and we should _not_ use this flag

  when: 'MuTect2' in workflowSteps

  script:
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar ${referenceMap['gatkHome']}/GenomeAnalysisTK.jar \
  -T MuTect2 \
  -R $genomeFile \
  --cosmic $cosmic \
  --dbsnp $dbsnp \
  -I:normal $bamNormal \
  -I:tumor $bamTumor \
  -U ALLOW_SEQ_DICT_INCOMPATIBILITY \
  -L \"$genInt\" \
  -XL hs37d5 \
  -XL NC_007605 \
  -o ${gen_int}_${idSampleNormal}_${idSampleTumor}.vcf
  """
}

if ('MuTect2' in workflowSteps) {
  if (verbose) {mutect2Output = mutect2Output.view {"MuTect2 output: $it"}}
  mutect2VCF = mutect2Output.groupTuple(by:[0,1,2,3,4])
  if (verbose) {mutect2VCF = mutect2VCF.view {"mutect2VCF: $it"}}
}

process RunFreeBayes {
  tag {idSampleTumor + "-" + gen_int}

  input:
    set idPatient, gender, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), genInt, gen_int from bamsFFB
    file genomeFile from file(referenceMap['genomeFile'])

  output:
    set val("FreeBayes"), idPatient, gender, idSampleNormal, idSampleTumor, val("${gen_int}_${idSampleNormal}_${idSampleTumor}"), file("${gen_int}_${idSampleNormal}_${idSampleTumor}.vcf") into freebayesOutput

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
    $bamNormal > ${gen_int}_${idSampleNormal}_${idSampleTumor}.vcf
  """
}

if ('FreeBayes' in workflowSteps) {
  if (verbose) {freebayesOutput = freebayesOutput.view {"FreeBayes output: $it"}}
  freebayesVCF = freebayesOutput.groupTuple(by:[0,1,2,3,4])
  if (verbose) {freebayesVCF = freebayesVCF.view {"freebayesVCF: $it"}}
}

process RunVardict {
  tag {idSampleTumor + "-" + gen_int}

  input:
    set idPatient, gender, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), genInt, gen_int from bamsFVD
    file genomeFile from file(referenceMap['genomeFile'])
    file genomeIndex from file(referenceMap['genomeIndex'])
    file genomeDict from file(referenceMap['genomeDict'])

  output:
    set val("VarDict"), idPatient, gender, idSampleNormal, idSampleTumor, val("${gen_int}_${idSampleNormal}_${idSampleTumor}"), file("${gen_int}_${idSampleNormal}_${idSampleTumor}.out") into vardictOutput

  when: 'VarDict' in workflowSteps

  script:
  """
  ${referenceMap['vardictHome']}/vardict.pl \
  -G $genomeFile \
  -f 0.01 -N $bamTumor \
  -b "$bamTumor|$bamNormal" \
  -z 1 -F 0x500 \
  -c 1 -S 2 -E 3 -g 4 \
  -R $genInt > ${gen_int}_${idSampleNormal}_${idSampleTumor}.out
  """
}

if ('VarDict' in workflowSteps) {
  if (verbose) {vardictOutput = vardictOutput.view {"VarDict output: $it"}}
  vardictVCF = vardictOutput.groupTuple(by:[0,1,2,3,4])
  if (verbose) {vardictVCF = vardictVCF.view {"vardictVCF: $it"}}
}

// we are merging the VCFs that are called separatelly for different intervals
// so we can have a single sorted VCF containing all the calls for a given caller
if ('HaplotypeCaller' in workflowSteps || 'MuTect1' in workflowSteps || 'MuTect2' in workflowSteps || 'FreeBayes' in workflowSteps || 'VarDict' in workflowSteps) {
  vcfsToMerge = hcVCF.mix(mutect1VCF, mutect2VCF, freebayesVCF, vardictVCF)
  if (verbose) {vcfsToMerge = vcfsToMerge.view {"VCFs To be merged: $it"}}
} else {
  vcfsToMerge.close()
}

process ConcatVCF {
  tag {variantCaller == 'HaplotypeCaller' ? variantCaller + "-" + idSampleNormal : variantCaller + "-" + idSampleNormal + "-" + idSampleTumor}

  publishDir "${directoryMap["$variantCaller"]}", mode: 'copy'

  input:
    set variantCaller, idPatient, gender, idSampleNormal, idSampleTumor, tag, file(vcFiles) from vcfsToMerge

  output:
    set variantCaller, idPatient, gender, idSampleNormal, idSampleTumor, file("*.vcf") into vcfConcatenated

  when: 'HaplotypeCaller' in workflowSteps || 'MuTect1' in workflowSteps || 'MuTect2' in workflowSteps || 'FreeBayes' in workflowSteps || 'VarDict' in workflowSteps

  script:
  vcfFiles = vcFiles.collect{"-V $it"}.join(' ')
  outputFile = (variantCaller == 'HaplotypeCaller' ? "${variantCaller}_${idPatient}_${idSampleNormal}.vcf" : "${variantCaller}_${idPatient}_${idSampleNormal}_${idSampleTumor}.vcf")

  if (variantCaller == 'VarDict')
    """
    for i in $vcFiles ;do
      cat \$i | ${referenceMap['vardictHome']}/VarDict/testsomatic.R >> testsomatic.out
    done
    ${referenceMap['vardictHome']}/VarDict/var2vcf_somatic.pl \
    -f 0.01 \
    -N "${idPatient}_${idSampleNormal}_${idSampleTumor}" testsomatic.out > $outputFile
    """

  else
    """
    java -Xmx${task.memory.toGiga()}g \
    -cp ${referenceMap['gatkHome']}/GenomeAnalysisTK.jar \
    org.broadinstitute.gatk.tools.CatVariants \
    --reference ${referenceMap['genomeFile']} \
    $vcfFiles \
    --outputFile $outputFile
    """
}

if ('HaplotypeCaller' in workflowSteps || 'MuTect1' in workflowSteps || 'MuTect2' in workflowSteps || 'FreeBayes' in workflowSteps || 'VarDict' in workflowSteps) {
   if (verbose) {vcfConcatenated = vcfConcatenated.view {"VCF concatenated: $it"}}
}

process RunStrelka {
  tag {idSampleTumor}

  publishDir directoryMap['Strelka'], mode: 'copy'

  input:
    set idPatient, gender, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from bamsForStrelka
    file genomeFile from file(referenceMap['genomeFile'])
    file genomeIndex from file(referenceMap['genomeIndex'])
    file genomeDict from file(referenceMap['genomeDict'])
    file strelkaCFG from file(referenceMap['strelkaCFG'])

  output:
    set val("Strelka"), idPatient, gender, idSampleNormal, idSampleTumor, file("strelka/results/*.vcf") into strelkaOutput

  when: 'Strelka' in workflowSteps

  script:
  """
  tumorPath=`readlink $bamTumor`
  normalPath=`readlink $bamNormal`
  ${referenceMap['strelkaHome']}/bin/configureStrelkaWorkflow.pl \
  --tumor \$tumorPath \
  --normal \$normalPath \
  --ref $genomeFile \
  --config $strelkaCFG \
  --output-dir strelka

  cd strelka

  make -j $task.cpus
  """
}

if ('Strelka' in workflowSteps) {
  if (verbose) {strelkaOutput = strelkaOutput.view {"Strelka output: $it"}}
}

process RunManta {
  tag {idSampleTumor}

  publishDir directoryMap['Manta'], mode: 'copy'

  input:
    set idPatient, gender, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from bamsForManta
    file mantaRef from file(referenceMap['mantaRef'])
    file mantaIndex from file(referenceMap['mantaIndex'])

  output:
    set val("Manta"), idPatient, gender, idSampleNormal, idSampleTumor, file("${idSampleNormal}_${idSampleTumor}.somaticSV.vcf"),file("${idSampleNormal}_${idSampleTumor}.candidateSV.vcf"),file("${idSampleNormal}_${idSampleTumor}.diploidSV.vcf"),file("${idSampleNormal}_${idSampleTumor}.candidateSmallIndels.vcf") into mantaOutput

  when: 'Manta' in workflowSteps

  //NOTE: Manta is very picky about naming and reference indexes, the input bam should not contain too many _ and the reference index must be generated using a supported samtools version.
  //Moreover, the bam index must be named .bam.bai, otherwise it will not be recognized
  script:
  """
  set -eo pipefail
  samtools view -H $bamNormal | grep -v hs37d5 | grep -v NC_007605 | samtools reheader - $bamNormal > Normal.bam
  samtools index Normal.bam

  samtools view -H $bamTumor | grep -v hs37d5 | grep -v NC_007605 | samtools reheader - $bamTumor > Tumor.bam
  samtools index Tumor.bam

  configManta.py --normalBam Normal.bam --tumorBam Tumor.bam --reference $mantaRef --runDir MantaDir
  python MantaDir/runWorkflow.py -m local -j $task.cpus
  gunzip -c MantaDir/results/variants/somaticSV.vcf.gz > ${idSampleNormal}_${idSampleTumor}.somaticSV.vcf
  gunzip -c MantaDir/results/variants/candidateSV.vcf.gz > ${idSampleNormal}_${idSampleTumor}.candidateSV.vcf
  gunzip -c MantaDir/results/variants/diploidSV.vcf.gz > ${idSampleNormal}_${idSampleTumor}.diploidSV.vcf
  gunzip -c MantaDir/results/variants/candidateSmallIndels.vcf.gz > ${idSampleNormal}_${idSampleTumor}.candidateSmallIndels.vcf
  """
}

if ('Manta' in workflowSteps) {
  if (verbose) {mantaOutput = mantaOutput.view {"Manta output: $it"}}
}

// Run commands and code from Malin Larsson
// Based on Jesper Eisfeldt's code
process RunAlleleCount {
  tag {idSampleTumor}

  input:
    set idPatient, gender, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from bamsForAscat
    file genomeFile from file(referenceMap['genomeFile'])
    file genomeIndex from file(referenceMap['genomeIndex'])
    file genomeDict from file(referenceMap['genomeDict'])
    file acLoci from file(referenceMap['acLoci'])

  output:
    set idPatient, gender, idSampleNormal, idSampleTumor, file("${idSampleNormal}.alleleCount"), file("${idSampleTumor}.alleleCount") into alleleCountOutput

  when: 'Ascat' in workflowSteps

  script:
  """
  alleleCounter -l $acLoci -r $genomeFile -b $bamNormal -o ${idSampleNormal}.alleleCount;
  alleleCounter -l $acLoci -r $genomeFile -b $bamTumor -o ${idSampleTumor}.alleleCount;
  """
}

// R script from Malin Larssons bitbucket repo:
// https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
process RunConvertAlleleCounts {
  tag {idSampleTumor}

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
  tag {idSampleTumor}

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

if ('Ascat' in workflowSteps) {
  if (verbose) {ascatOutput = ascatOutput.view {"Ascat output: $it"}}
}

// reportsForMultiQC = Channel.create()
//
// if ('MultiQC' in workflowSteps) {
//   reportsForMultiQC = reportsForMultiQC.mix(fastQCreport).flatten().toList()
//   if (verbose) {reportsForMultiQC = reportsForMultiQC.view {"Reports for MultiQC: $it"}}
// }
//
// process RunMultiQC {
//   tag {"MultiQC"}
//
//   publishDir directoryMap['MultiQC'], mode: 'copy'
//
//   input:
//     file ('*') from reportsForMultiQC
//
//   output:
//     set file("*multiqc_report.html"), file("*multiqc_data") into multiQCReport
//
//     when: 'MultiQC' in workflowStepsworkflowSteps
//
//   script:
//   """
//   multiqc -f -v .
//   """
// }
// if ('MultiQC' in workflowSteps) {
//   if (verbose) {multiQCReport = multiQCReport.view {"MultiQC report: $it"}}
// }

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
    'genomeFile'  : params.genome,      // genome reference
    'genomeIndex' : params.genomeIndex, // genome reference index
    'genomeDict'  : params.genomeDict,  // genome reference dictionary
    'kgIndels'    : params.kgIndels,    // 1000 Genomes SNPs
    'kgIndex'     : params.kgIndex,     // 1000 Genomes SNPs index
    'dbsnp'       : params.dbsnp,       // dbSNP
    'dbsnpIndex'  : params.dbsnpIndex,  // dbSNP index
    'millsIndels' : params.millsIndels, // Mill's Golden set of SNPs
    'millsIndex'  : params.millsIndex,  // Mill's Golden set index
    'cosmic'      : params.cosmic,      // cosmic vcf file with VCF4.1 header
    'cosmicIndex' : params.cosmicIndex, // cosmic vcf file index
    'intervals'   : params.intervals,   // intervals file for spread-and-gather processes (usually chromosome chunks at centromeres)
    'mantaRef'    : params.mantaRef,    // copy of the genome reference file
    'mantaIndex'  : params.mantaIndex,  // reference index indexed with samtools/0.1.19
    'acLoci'      : params.acLoci,      // loci file for ascat
    'picardHome'  : params.picardHome,  // path to Picard
    'gatkHome'    : params.gatkHome,    // path to Gatk
    'mutect1Home' : params.mutect1Home, // path to MuTect1
    'vardictHome' : params.vardictHome, // path to VarDict
    'strelkaHome' : params.strelkaHome, // path to Strelka
    'strelkaCFG'  : params.strelkaCFG   // Strelka Config file
  ]
}

def defineDirectoryMap() {
  return [
    'nonRealigned'    : 'Preprocessing/NonRealigned',
    'recalibrated'    : 'Preprocessing/Recalibrated',
    'MuTect1'         : 'VariantCalling/MuTect1',
    'MuTect2'         : 'VariantCalling/MuTect2',
    'FreeBayes'       : 'VariantCalling/FreeBayes',
    'VarDict'         : 'VariantCalling/VarDict',
    'Strelka'         : 'VariantCalling/Strelka',
    'HaplotypeCaller' : 'VariantCalling/HaplotypeCaller',
    'Manta'           : 'VariantCalling/Manta',
    'Ascat'           : 'VariantCalling/Ascat',
    'MultiQC'         : 'Reports/MultiQC'
  ]
}

def defineStepList() {
  return [
    'preprocessing',
    'realign',
    'skipPreprocessing',
    'MuTect1',
    'MuTect2',
    'FreeBayes',
    'VarDict',
    'Strelka',
    'HaplotypeCaller',
    'Manta',
    'Ascat',
    'MultiQC'
  ]
}

def checkReferenceMap(referenceMap) {
  referenceDefined = true
  referenceMap.each{ //Loop through all the references files to check their existence
    referenceFile, fileToCheck ->
    test = checkRefExistence(referenceFile, fileToCheck)
    !(test) ? referenceDefined = false : ""
  }
  return (referenceDefined ? true : false)
}

def checkStepList(stepsList, realStepsList) {
  stepCorrect = true
  stepsList.each{ // Loop through all the possible steps check their existence and spelling
    test = checkStepExistence(it, realStepsList)
    !(test) ? stepCorrect = false : ""
  }
  return (stepCorrect ? true : false)
}

def checkRefExistence(referenceFile, fileToCheck) { // Check file existence
  try {assert file(fileToCheck).exists()}
  catch (AssertionError ae) {
    log.info  "Missing references: $referenceFile $fileToCheck"
    return false
  }
  return true
}

def checkStepExistence(step, list) { // Check step existence
  try {assert list.contains(step)}
  catch (AssertionError ae) {
    println("Unknown parameter: $step")
    return false
  }
  return true
}

def checkFileExistence(it) { // Check file existence
  try {assert file(it).exists()}
  catch (AssertionError ae) {
    exit 1, "Missing file in TSV file: $it, see --help for more information"
  }
}

def extractFastqFiles(tsvFile) { // Channeling the TSV file containing FASTQ. Format is: "subject gender status sample lane fastq1 fastq2"
  fastqFiles = Channel
    .from(tsvFile.readLines())
    .map{line ->
      list       = line.split()
      idPatient  = list[0]
      gender     = list[1]
      status     = list[2]
      idSample   = list[3]
      idRun      = list[4]
      temp1      = list[5]
      temp2      = list[6]

      // When testing workflow from github, paths to FASTQ files start from workflow.projectDir and not workflow.launchDir
      if ((workflow.commitId) && (params.test)) {
        fastqFile1 = file("$workflow.projectDir/$temp1")
        fastqFile2 = file("$workflow.projectDir/$temp2")
      } else {
        fastqFile1 = file("$temp1")
        fastqFile2 = file("$temp2")
      }

      checkFileExistence(fastqFile1)
      checkFileExistence(fastqFile2)

      [idPatient, gender, status, idSample, idRun, fastqFile1, fastqFile2]
    }
  return fastqFiles
}

def extractBamFiles(tsvFile) { // Channeling the TSV file containing BAM. Format is: "subject gender status sample bam bai"
  bamFiles = Channel
    .from(tsvFile.readLines())
    .map{line ->
      list      = line.split()
      idPatient = list[0]
      gender    = list[1]
      status    = list[2]
      idSample  = list[3]
      bamFile   = file(list[4])
      baiFile   = file(list[5])

      checkFileExistence(bamFile)
      checkFileExistence(baiFile)

      [ idPatient, gender, status, idSample, bamFile, baiFile ]
    }
  return bamFiles
}

def retreiveStatus(bamChannel) {
  bamChannel = bamChannel.map {
    idPatient, gender, tag, bam, bai ->
    array = tag.split("_")
    status   = array[1]
    idSample = array[0]
    [idPatient, gender, status, idSample, bam, bai]
  }
  return bamChannel
}

def generateIntervalsForVC(bams, gI) {
  bamsForVC = Channel.create()
  vcIntervals = Channel.create()
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
  if (workflow.commitId) {
    log.info "Git info    : $workflow.repository - $workflow.revision [$workflow.commitId]"
  } else {
    log.info "  revision  : $revision"
  }
}

workflow.onComplete { // Display complete message
  log.info "N E X T F L O W ~ $workflow.nextflow.version - $workflow.nextflow.build"
  log.info "CANCER ANALYSIS WORKFLOW ~ $version - revision: $revision"
  log.info "Command Line: $workflow.commandLine"
  log.info "Project Dir : $workflow.projectDir"
  log.info "Launch Dir  : $workflow.launchDir"
  log.info "Work Dir    : $workflow.workDir"
  log.info "Steps       : " + workflowSteps.join(", ")
  log.info "Completed at: $workflow.complete"
  log.info "Duration    : $workflow.duration"
  log.info "Success     : $workflow.success"
  log.info "Exit status : $workflow.exitStatus"
  log.info "Error report: " + (workflow.errorReport ?: '-')
}

workflow.onError { // Display error message
  log.info "N E X T F L O W ~ version $workflow.nextflow.version [$workflow.nextflow.build]"
  if (workflow.commitId) {
    log.info "CANCER ANALYSIS WORKFLOW ~ $version - $workflow.revision [$workflow.commitId]"
  } else {
    log.info "CANCER ANALYSIS WORKFLOW ~ $version - revision: $revision"
  }
  log.info "Workflow execution stopped with the following message: " + workflow.errorMessage
}
