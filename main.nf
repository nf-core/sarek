#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
kate: syntax groovy; space-indent on; indent-width 4;
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
 - MapReads - Map reads
 - MergeBams - Merge BAMs if multilane samples
 - MarkDuplicates - Mark Duplicates
 - RealignerTargetCreator - Create realignment target intervals
 - IndelRealigner - Realign BAMs as T/N pair
 - CreateRecalibrationTable - Create Recalibration Table
 - RecalibrateBam - Recalibrate Bam
 - RunSamtoolsStats - Run Samtools stats on recalibrated BAM files
 - RunHaplotypecaller - Run HaplotypeCaller for GermLine Variant Calling (Parallelized processes)
 - RunMutect1 - Run MuTect1 for Variant Calling (Parallelized processes)
 - RunMutect2 - Run MuTect2 for Variant Calling (Parallelized processes)
 - RunFreeBayes - Run FreeBayes for Variant Calling (Parallelized processes)
 - RunVardict - Run VarDict for Variant Calling (Parallelized processes)
 - ConcatVCF - Merge results from HaplotypeCaller, MuTect1, MuTect2 and VarDict
 - RunStrelka - Run Strelka for Variant Calling
 - RunManta - Run Manta for Structural Variant Calling
 - RunAlleleCount - Run AlleleCount to prepare for ASCAT
 - RunConvertAlleleCounts - Run convertAlleleCounts to prepare for ASCAT
 - RunAscat - Run ASCAT for CNV
 - RunSnpeff - Run snpEff for annotation of vcf files
 - RunVEP - Run VEP for annotation of vcf files
 - RunBcftoolsStats - Run BCFTools stats on vcf files
 - GenerateMultiQCconfig - Generate a config file for MultiQC
 - RunMultiQC - Run MultiQC for report and QC
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

version = '1.1'

if (!isAllowedParams(params)) {exit 1, "params is unknown, see --help for more information"}

if (params.help) {
  helpMessage(version, grabRevision())
  exit 1
}

if (params.version) {
  versionMessage(version, grabRevision())
  exit 1
}

if (!checkUppmaxProject()) {exit 1, 'No UPPMAX project ID found! Use --project <UPPMAX Project ID>'}

step = params.step.toLowerCase()
tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase()} : []

directoryMap = defineDirectoryMap()
referenceMap = defineReferenceMap()
stepList = defineStepList()
toolList = defineToolList()
verbose = params.verbose

if (!checkExactlyOne([params.test, params.sample, params.sampleDir]))
  exit 1, 'Please define which samples to work on by providing exactly one of the --test, --sample or --sampleDir options'
if (!checkReferenceMap(referenceMap)) {exit 1, 'Missing Reference file(s), see --help for more information'}
if (!checkParameterExistence(step, stepList)) {exit 1, 'Unknown step, see --help for more information'}
if (step.contains(',')) {exit 1, 'You can choose only one step, see --help for more information'}
if (!checkParameterList(tools,toolList)) {exit 1, 'Unknown tool(s), see --help for more information'}

tsvPath = ''
if (params.test) {
  if (params.genome == "GRCh37") {
    referenceMap.intervals = file("$workflow.projectDir/repeats/tiny.list")
  }
  testTsvPaths = [
    'preprocessing': "$workflow.projectDir/data/tsv/tiny.tsv",
    'realign': "$workflow.launchDir/$directoryMap.nonRealigned/nonRealigned.tsv",
    'recalibrate': "$workflow.launchDir/$directoryMap.nonRecalibrated/nonRecalibrated.tsv",
    'skippreprocessing': "$workflow.launchDir/$directoryMap.recalibrated/recalibrated.tsv"
  ]
  tsvPath = testTsvPaths[step]
} else if (params.sample) {
  tsvPath = params.sample
}

// Set up the fastqFiles and bamFiles channels. One of them remains empty
fastqFiles = Channel.empty()
bamFiles = Channel.empty()
if (tsvPath) {
  tsvFile = file(tsvPath)
  switch (step) {
    case 'preprocessing': fastqFiles = extractFastq(tsvFile); break
    case 'realign': bamFiles = extractBams(tsvFile); break
    case 'recalibrate': bamFiles = extractRecal(tsvFile); break
    case 'skippreprocessing': bamFiles = extractBams(tsvFile); break
    default: exit 1, "Unknown step $step"
  }
} else if (params.sampleDir) {
  if (step != 'preprocessing') exit 1, '--sampleDir does not support steps other than "preprocessing"'
  fastqFiles = extractFastqFromDir(params.sampleDir)
  tsvFile = params.sampleDir  // used in the reports
}
verbose ? fastqFiles = fastqFiles.view {"FASTQ files to preprocess: $it"} : ''
verbose ? bamFiles = bamFiles.view {"BAM files to process: $it"} : ''
startMessage(version, grabRevision())

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

(fastqFiles, fastqFilesforFastQC) = fastqFiles.into(2)

verbose ? fastqFilesforFastQC = fastqFilesforFastQC.view {"FASTQ files for FastQC: $it"} : ''

process RunFastQC {
  tag {idPatient + "-" + idRun}

  publishDir directoryMap.fastQC, mode: 'copy'

  input:
    set idPatient, gender, status, idSample, idRun, file(fastqFile1), file(fastqFile2) from fastqFilesforFastQC

  output:
    file "*_fastqc.{zip,html}" into fastQCreport

  when: step == 'preprocessing' && 'multiqc' in tools

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
    set file(genomeFile), file(bwaIndex) from Channel.value([referenceMap.genomeFile, referenceMap.bwaIndex])

  output:
    set idPatient, gender, status, idSample, idRun, file("${idRun}.bam") into mappedBam

  when: step == 'preprocessing'

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

verbose ? mappedBam = mappedBam.view {"BAM file to sort into group or single: $it"} : ''

// Sort bam whether they are standalone or should be merged
// Borrowed code from https://github.com/guigolab/chip-nf

singleBam = Channel.create()
groupedBam = Channel.create()
mappedBam.groupTuple(by:[0,1,2,3])
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

  when: step == 'preprocessing'

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

  publishDir '.', saveAs: { it == "${bam}.metrics" ? "$directoryMap.markDuplicatesQC/$it" : "$directoryMap.nonRealigned/$it" }, mode: 'copy'

  input:
    set idPatient, gender, status, idSample, file(bam) from mergedBam

  output:
    set idPatient, gender, val("${idSample}_${status}"), file("${idSample}_${status}.md.bam"), file("${idSample}_${status}.md.bai") into duplicates
    set idPatient, gender, status, idSample, val("${idSample}_${status}.md.bam"), val("${idSample}_${status}.md.bai") into markDuplicatesTSV
    file ("${bam}.metrics") into markDuplicatesReport

  when: step == 'preprocessing'

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
  "$idPatient\t$gender\t$status\t$idSample\t$directoryMap.nonRealigned/$bam\t$directoryMap.nonRealigned/$bai\n"
}.collectFile(
  name: 'nonRealigned.tsv', sort: true, storeDir: directoryMap.nonRealigned
)

// Create intervals for realignement using both tumor+normal as input
// Group the marked duplicates BAMs for intervals and realign by idPatient
// Grouping also by gender, to make a nicer channel
duplicatesGrouped = step == 'preprocessing' ? duplicates.groupTuple(by:[0,1]) : Channel.empty()

duplicatesGrouped = step == 'realign' ? bamFiles.map{
  idPatient, gender, status, idSample, bam, bai ->
  [idPatient, gender, "${idSample}_${status}", bam, bai]
}.groupTuple(by:[0,1]) : duplicatesGrouped

// The duplicatesGrouped channel is duplicated
// one copy goes to the RealignerTargetCreator process
// and the other to the IndelRealigner process
(duplicatesInterval, duplicatesRealign) = duplicatesGrouped.into(2)

verbose ? duplicatesInterval = duplicatesInterval.view {"BAMs for RealignerTargetCreator: $it"} : ''
verbose ? duplicatesRealign = duplicatesRealign.view {"BAMs to phase: $it"} : ''
verbose ? markDuplicatesReport = markDuplicatesReport.view {"MarkDuplicates report: $it"} : ''

// VCF indexes are added so they will be linked, and not re-created on the fly
//  -L "1:131941-141339" \

process RealignerTargetCreator {
  tag {idPatient}

  input:
    set idPatient, gender, idSample_status, file(bam), file(bai) from duplicatesInterval
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(knownIndels), file(knownIndelsIndex), file(intervals) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.knownIndels,
      referenceMap.knownIndelsIndex,
      referenceMap.intervals
    ])

  output:
    set idPatient, gender, file("${idPatient}.intervals") into intervals

  when: step == 'preprocessing' || step == 'realign'

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

verbose ? bamsAndIntervals = bamsAndIntervals.view {"BAMs and Intervals phased for IndelRealigner: $it"} : ''

// use nWayOut to split into T/N pair again
process IndelRealigner {
  tag {idPatient}

  input:
    set idPatient, gender, idSample_status, file(bam), file(bai), file(intervals) from bamsAndIntervals
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(knownIndels), file(knownIndelsIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.knownIndels,
      referenceMap.knownIndelsIndex])

  output:
    set idPatient, gender, file("*.real.bam"), file("*.real.bai") into realignedBam mode flatten

  when: step == 'preprocessing' || step == 'realign'

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

realignedBam = retrieveStatus(realignedBam)

verbose ? realignedBam = realignedBam.view {"Realigned BAM to CreateRecalibrationTable: $it"} : ''

process CreateRecalibrationTable {
  tag {idPatient + "-" + idSample}

  publishDir directoryMap.nonRecalibrated, mode: 'copy'

  input:
    set idPatient, gender, status, idSample, file(bam), file(bai) from realignedBam
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
    set idPatient, gender, status, idSample, file(bam), file(bai), file("${idSample}.recal.table") into recalibrationTable
    set idPatient, gender, status, idSample, val("${idSample}_${status}.md.real.bam"), val("${idSample}_${status}.md.real.bai"), val("${idSample}.recal.table") into recalibrationTableTSV

  when: step == 'preprocessing' || step == 'realign'

  script:
  known = knownIndels.collect{ "-knownSites $it" }.join(' ')
  """
  java -Xmx${task.memory.toGiga()}g \
  -Djava.io.tmpdir="/tmp" \
  -jar \$GATK_HOME/GenomeAnalysisTK.jar \
  -T BaseRecalibrator \
  -R $genomeFile \
  -I $bam \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  -knownSites $dbsnp \
  $known \
  -nct $task.cpus \
  -L $intervals \
  -l INFO \
  -o ${idSample}.recal.table
  """
}
// Creating a TSV file to restart from this step
recalibrationTableTSV.map { idPatient, gender, status, idSample, bam, bai, recalTable ->
  "$idPatient\t$gender\t$status\t$idSample\t$directoryMap.nonRecalibrated/$bam\t$directoryMap.nonRecalibrated/$bai\t\t$directoryMap.nonRecalibrated/$recalTable\n"
}.collectFile(
  name: 'nonRecalibrated.tsv', sort: true, storeDir: directoryMap.nonRecalibrated
)

recalibrationTable = step == 'recalibrate' ? bamFiles : recalibrationTable

verbose ? recalibrationTable = recalibrationTable.view {"Base recalibrated table for RecalibrateBam: $it"} : ''

process RecalibrateBam {
  tag {idPatient + "-" + idSample}

  publishDir directoryMap.recalibrated, mode: 'copy'

  input:
    set idPatient, gender, status, idSample, file(bam), file(bai), recalibrationReport from recalibrationTable
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(intervals) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.intervals,
    ])

  output:
    set idPatient, gender, status, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bai") into recalibratedBam
    set idPatient, gender, status, idSample, val("${idSample}.recal.bam"), val("${idSample}.recal.bai") into recalibratedBamTSV

  when: step != 'skippreprocessing'

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
recalibratedBamTSV.map { idPatient, gender, status, idSample, bam, bai ->
  "$idPatient\t$gender\t$status\t$idSample\t$directoryMap.recalibrated/$bam\t$directoryMap.recalibrated/$bai\n"
}.collectFile(
  name: 'recalibrated.tsv', sort: true, storeDir: directoryMap.recalibrated
)

recalibratedBam = step == 'skippreprocessing' ? bamFiles : recalibratedBam

verbose ? recalibratedBam = recalibratedBam.view {"Recalibrated BAM for variant Calling: $it"} : ''

(recalibratedBam, recalibratedBamForStats) = recalibratedBam.into(2)

process RunSamtoolsStats {
  tag {idPatient + "-" + idSample}

  publishDir directoryMap.samtoolsStats, mode: 'copy'

  input:
    set idPatient, gender, status, idSample, file(bam), file(bai) from recalibratedBamForStats

  output:
    file ("${bam}.samtools.stats.out") into recalibratedBamReport

    when: 'multiqc' in tools

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
  .choice(bamsTumor, bamsNormal) {it[2] == 0 ? 1 : 0}

// Removing status because not relevant anymore
bamsNormal = bamsNormal.map { idPatient, gender, status, idSample, bam, bai -> [idPatient, gender, idSample, bam, bai] }
verbose ? bamsNormal = bamsNormal.view {"Normal BAM for variant Calling: $it"} : ''

bamsTumor = bamsTumor.map { idPatient, gender, status, idSample, bam, bai -> [idPatient, gender, idSample, bam, bai] }
verbose ? bamsTumor = bamsTumor.view {"Tumor BAM for variant Calling: $it"} : ''

// We know that MuTect2 (and other somatic callers) are notoriously slow. To speed them up we are chopping the reference into
// smaller pieces at centromeres (see repeats/centromeres.list), do variant calling by this intervals, and re-merge the VCFs.
// Since we are on a cluster, this can parallelize the variant call processes, and push down the variant call wall clock time significanlty.

// in fact we need two channels: one for the actual genomic region, and an other for names
// without ":", as nextflow is not happy with them (will report as a failed process).
// For region 1:1-2000 the output file name will be something like 1_1-2000_Sample_name.mutect2.vcf
// from the "1:1-2000" string make ["1:1-2000","1_1-2000"]

// define intervals file by --intervals
intervals = Channel.from(file(referenceMap.intervals).readLines())
gI = intervals.map{[it,it.replaceFirst(/\:/,'_')]}

(bamsNormalTemp, bamsNormal, gI) = generateIntervalsForVC(bamsNormal, gI)
(bamsTumorTemp, bamsTumor, gI) = generateIntervalsForVC(bamsTumor, gI)

// HaplotypeCaller
bamsFHC = bamsNormalTemp.mix(bamsTumorTemp)
verbose ? bamsFHC = bamsFHC.view {"Bams with Intervals for HaplotypeCaller: $it"} : ''

(bamsNormalTemp, bamsNormal) = bamsNormal.into(2)
(bamsTumorTemp, bamsTumor) = bamsTumor.into(2)

bamsNormalTemp = bamsNormalTemp.map { idPatient, gender, idSample, bam, bai -> [idPatient, gender, 0, idSample, bam, bai] }
bamsTumorTemp = bamsTumorTemp.map { idPatient, gender, idSample, bam, bai -> [idPatient, gender, 1, idSample, bam, bai] }

bamsForAscat = Channel.create()
bamsForAscat = bamsNormalTemp.mix(bamsTumorTemp)
verbose ? bamsForAscat = bamsForAscat.view {"Bams for Ascat: $it"} : ''

bamsAll = bamsNormal.spread(bamsTumor)
// Since idPatientNormal and idPatientTumor are the same
// It's removed from bamsAll Channel (same for genderNormal)
// /!\ It is assumed that every sample are from the same patient
bamsAll = bamsAll.map {
  idPatientNormal, genderNormal, idSampleNormal, bamNormal, baiNormal, idPatientTumor, genderTumor, idSampleTumor, bamTumor, baiTumor ->
  [idPatientNormal, genderNormal, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor]
}
verbose ? bamsAll = bamsAll.view {"Mapped Recalibrated BAM for variant Calling: $it"} : ''

// MuTect1
(bamsFMT1, bamsAll, gI) = generateIntervalsForVC(bamsAll, gI)
verbose ? bamsFMT1 = bamsFMT1.view {"Bams with Intervals for MuTect1: $it"} : ''

// MuTect2
(bamsFMT2, bamsAll, gI) = generateIntervalsForVC(bamsAll, gI)
verbose ? bamsFMT2 = bamsFMT2.view {"Bams with Intervals for MuTect2: $it"} : ''

// FreeBayes
(bamsFFB, bamsAll, gI) = generateIntervalsForVC(bamsAll, gI)
verbose ? bamsFFB = bamsFFB.view {"Bams with Intervals for FreeBayes: $it"} : ''

// VarDict
(bamsFVD, bamsAll, gI) = generateIntervalsForVC(bamsAll, gI)
verbose ? bamsFVD = bamsFVD.view {"Bams with Intervals for VarDict: $it"} : ''

(bamsForManta, bamsForStrelka) = bamsAll.into(2)

verbose ? bamsForManta = bamsForManta.view {"Bams for Manta: $it"} : ''

verbose ? bamsForStrelka = bamsForStrelka.view {"Bams for Strelka: $it"} : ''

process RunHaplotypecaller {
  tag {idPatient + "-" + idSample + "-" + gen_int}

  input:
    set idPatient, gender, idSample, file(bam), file(bai), genInt, gen_int from bamsFHC //Are these values `ped to bamNormal already?
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(dbsnp), file(dbsnpIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.dbsnp,
      referenceMap.dbsnpIndex
    ])

  output:
    set val("haplotypecaller"), idPatient, gender, idSample, val("${gen_int}_${idSample}"), file("${gen_int}_${idSample}.g.vcf") into hcVCF

  when: 'haplotypecaller' in tools

  script:
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar \$GATK_HOME/GenomeAnalysisTK.jar \
  -T HaplotypeCaller \
  --emitRefConfidence GVCF \
  -pairHMM LOGLESS_CACHING \
  -R $genomeFile \
  --dbsnp $dbsnp \
  -I $bam \
  -L \"$genInt\" \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  -o ${gen_int}_${idSample}.g.vcf
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
    set val("mutect1"), idPatient, gender, idSampleNormal, idSampleTumor, val("${gen_int}_${idSampleTumor}_vs_${idSampleNormal}"), file("${gen_int}_${idSampleTumor}_vs_${idSampleNormal}.vcf") into mutect1Output

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
  -L \"$genInt\" \
  --disable_auto_index_creation_and_locking_when_reading_rods \
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
    set val("mutect2"), idPatient, gender, idSampleNormal, idSampleTumor, val("${gen_int}_${idSampleTumor}_vs_${idSampleNormal}"), file("${gen_int}_${idSampleTumor}_vs_${idSampleNormal}.vcf") into mutect2Output

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
  -L \"$genInt\" \
  -o ${gen_int}_${idSampleTumor}_vs_${idSampleNormal}.vcf
  """
}

mutect2Output = mutect2Output.groupTuple(by:[0,1,2,3,4])
verbose ? mutect2Output = mutect2Output.view {"MuTect2 output: $it"} : ''

process RunFreeBayes {
  tag {idPatient + "-" + idSampleTumor + "-" + gen_int}

  input:
    set idPatient, gender, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), genInt, gen_int from bamsFFB
    file(genomeFile) from Channel.value(referenceMap.genomeFile)

  output:
    set val("freebayes"), idPatient, gender, idSampleNormal, idSampleTumor, val("${gen_int}_${idSampleTumor}_vs_${idSampleNormal}"), file("${gen_int}_${idSampleTumor}_vs_${idSampleNormal}.vcf") into freebayesOutput

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
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set val("vardict"), idPatient, gender, idSampleNormal, idSampleTumor, val("${gen_int}_${idSampleTumor}_vs_${idSampleNormal}"), file("${gen_int}_${idSampleTumor}_vs_${idSampleNormal}.out") into vardictOutput

  when: 'vardict' in tools

  script:
  """
  vardict.pl \
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
  tag {variantCaller == 'haplotypecaller' ? idPatient + "-" + variantCaller + "-" + idSampleNormal : idPatient + "-" + variantCaller + "-" + idSampleNormal + "-" + idSampleTumor}

  publishDir "${directoryMap."$variantCaller"}", mode: 'copy'

  input:
    set variantCaller, idPatient, gender, idSampleNormal, idSampleTumor, tag, file(vcFiles) from vcfsToMerge
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set variantCaller, idPatient, gender, idSampleNormal, idSampleTumor, file("*.vcf.gz") into vcfConcatenated

  when: 'haplotypecaller' in tools || 'mutect1' in tools || 'mutect2' in tools || 'freebayes' in tools || 'vardict' in tools

  script:
  outputFile = variantCaller == 'haplotypecaller' ? "${variantCaller}_${idSampleNormal}.vcf" : "${variantCaller}_${idSampleTumor}_vs_${idSampleNormal}.vcf"
  vcfFiles = vcFiles.collect{" $it"}.join(' ')

  if (variantCaller == 'vardict')
    """
    for i in $vcFiles ;do
      cat \$i | testsomatic.R >> testsomatic.out
    done
    var2vcf_paired.pl \
    -f 0.01 \
    -N "${idSampleTumor}_vs_${idSampleNormal}" testsomatic.out > $outputFile
    gzip -v $outputFile
    """

  else if (variantCaller == 'mutect2' || variantCaller == 'mutect1' || variantCaller == 'haplotypecaller' || variantCaller == 'freebayes')
    """
    # first make a header from one of the VCF intervals
    # get rid of interval information only from the GATK command-line, but leave the rest
    awk '/^#/{print}' `ls *vcf| head -1` | \
    awk '!/GATKCommandLine/{print}/GATKCommandLine/{for(i=1;i<=NF;i++){if(\$i!~/intervals=/ && \$i !~ /out=/){printf("%s ",\$i)}}printf("\\n")}' \
    > header

    ## concatenate calls
    ( cat header && awk '!/^#/' $vcfFiles ) > unsorted.vcf
    mkfifo out.vcf
    gzip < out.vcf > ${outputFile}.gz &
    gzpid=\$!
    java -jar \${PICARD_HOME}/picard.jar SortVcf I=unsorted.vcf O=out.vcf SEQUENCE_DICTIONARY=$genomeDict
    wait \$gzpid
    """
}

verbose ? vcfConcatenated = vcfConcatenated.view {"VCF concatenated: $it"} : ''

process RunStrelka {
  tag {idPatient + "-" + idSampleTumor}

  publishDir directoryMap.strelka, mode: 'copy'

  input:
    set idPatient, gender, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from bamsForStrelka
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set val("strelka"), idPatient, gender, idSampleNormal, idSampleTumor, file("*.vcf") into strelkaOutput

  when: 'strelka' in tools

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

  publishDir directoryMap.manta, mode: 'copy'

  input:
    set idPatient, gender, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from bamsForManta
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex
    ])

  output:
    set val("manta"), idPatient, gender, idSampleNormal, idSampleTumor, file("Manta_${idSampleTumor}_vs_${idSampleNormal}.somaticSV.vcf"),file("Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSV.vcf"),file("Manta_${idSampleTumor}_vs_${idSampleNormal}.diploidSV.vcf"),file("Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSmallIndels.vcf") into mantaOutput

  when: 'manta' in tools

  script:
  """
  ln -s $bamNormal Normal.bam
  ln -s $bamTumor Tumor.bam
  ln -s $baiNormal Normal.bam.bai
  ln -s $baiTumor Tumor.bam.bai

  \$MANTA_INSTALL_PATH/bin/configManta.py --normalBam Normal.bam --tumorBam Tumor.bam --reference $genomeFile --runDir MantaDir
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
    set file(acLoci), file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.acLoci,
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set idPatient, gender, status, idSample, file("${idSample}.alleleCount") into alleleCountOutput

  when: 'ascat' in tools

  script:
  """
  alleleCounter -l $acLoci -r $genomeFile -b $bam -o ${idSample}.alleleCount;
  """
}

verbose ? alleleCountOutput = alleleCountOutput.view {"alleleCount output: $it"} : ''

alleleCountNormal = Channel.create()
alleleCountTumor = Channel.create()

alleleCountOutput
  .choice(alleleCountTumor, alleleCountNormal) {it[2] == 0 ? 1 : 0}

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

  publishDir directoryMap.ascat, mode: 'copy'

  input:
    set idPatient, gender, idSampleNormal, idSampleTumor, file(alleleCountNormal), file(alleleCountTumor) from alleleCountOutput

  output:
    set idPatient, gender, idSampleNormal, idSampleTumor, file("${idSampleNormal}.BAF"), file("${idSampleNormal}.LogR"), file("${idSampleTumor}.BAF"), file("${idSampleTumor}.LogR") into convertAlleleCountsOutput

  when: 'ascat' in tools

  script:
  """
  convertAlleleCounts.r $idSampleTumor $alleleCountTumor $idSampleNormal $alleleCountNormal $gender
  """
}

// R scripts from Malin Larssons bitbucket repo:
// https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
process RunAscat {
  tag {idPatient + "-" + idSampleTumor}

  publishDir directoryMap.ascat, mode: 'copy'

  input:
    set idPatient, gender, idSampleNormal, idSampleTumor, file(bafNormal), file(logrNormal), file(bafTumor), file(logrTumor) from convertAlleleCountsOutput

  output:
    set val("ascat"), idPatient, gender, idSampleNormal, idSampleTumor, file("${idSampleTumor}.tumour.png"), file("${idSampleTumor}.germline.png"), file("${idSampleTumor}.LogR.PCFed.txt"), file("${idSampleTumor}.BAF.PCFed.txt"), file("${idSampleTumor}.ASPCF.png"), file("${idSampleTumor}.ASCATprofile.png"), file("${idSampleTumor}.aberrationreliability.png"), file("${idSampleTumor}.rawprofile.png"), file("${idSampleTumor}.sunrise.png"), file("${idSampleTumor}.cnvs.txt"),file("${idSampleTumor}.purityploidy.txt") into ascatOutput

  when: 'ascat' in tools

  script:
  """
  run_ascat.r $bafTumor $logrTumor $bafNormal $logrNormal $idSampleTumor $baseDir
  """
}

verbose ? ascatOutput = ascatOutput.view {"Ascat output: $it"} : ''

// process MergeVCF {
//   tag {idPatient + "-" + idSample}
//
//   input:
//     set idPatient, gender, status, idSample, file(vcf) from ???
//
//   output:
//     set idPatient, gender, status, idSample, file ("???") into vcfMerged
//
//     'mutect1' in tools || 'mutect2' in tools
//
//     script:
//     """
//
//     """
// }

vcfMerged = Channel.create()
vcfNotMerged = Channel.create()

vcfConcatenated
  .choice(vcfMerged, vcfNotMerged) {it[0] == 'mutect1' ? 0 : 1}

(strelkaAllIndels, strelkaAllSNVS, strelkaPAssedIndels, strelkaPAssedSNVS) = strelkaOutput.into(4)

strelkaAllIndels = strelkaAllIndels
  .map {
    variantcaller, idPatient, gender, idSampleNormal, idSampleTumor, vcf ->
    [variantcaller, idPatient, gender, idSampleNormal, idSampleTumor, vcf[0]]
}

strelkaAllSNVS = strelkaAllSNVS
  .map {
    variantcaller, idPatient, gender, idSampleNormal, idSampleTumor, vcf ->
    [variantcaller, idPatient, gender, idSampleNormal, idSampleTumor, vcf[1]]
}

strelkaPAssedIndels = strelkaPAssedIndels
  .map {
    variantcaller, idPatient, gender, idSampleNormal, idSampleTumor, vcf ->
    [variantcaller, idPatient, gender, idSampleNormal, idSampleTumor, vcf[2]]
}

strelkaPAssedSNVS = strelkaPAssedSNVS
  .map {
    variantcaller, idPatient, gender, idSampleNormal, idSampleTumor, vcf ->
    [variantcaller, idPatient, gender, idSampleNormal, idSampleTumor, vcf[3]]
}

vcfMerged = vcfMerged.mix(strelkaAllIndels, strelkaAllSNVS, strelkaPAssedIndels, strelkaPAssedSNVS)

verbose ? vcfMerged = vcfMerged.view {"VCF for Annotation: $it"} : ''

vcfForSnpeff = Channel.create()
vcfForBCF = Channel.create()
vcfForVep = Channel.create()

(vcfForBCF, vcfMerged) = vcfMerged.into(2)

process RunBcftoolsStats {
  tag {variantCaller + "-" + idPatient + "-" + idSample}

  publishDir directoryMap.bcftoolsStats, mode: 'copy'

  input:
    set variantCaller, idPatient, gender, idSampleNormal, idSampleTumor, file(vcf) from vcfForBCF

  output:
    file ("${vcf.baseName}.bcf.tools.stats.out") into bcfReport

  when: 'multiqc' in tools

  script:
  """
  bcftools stats $vcf > ${vcf.baseName}.bcf.tools.stats.out
  """
}

(vcfForSnpeff, vcfForVep) = vcfMerged.into(2)

process RunSnpeff {
  tag {variantCaller + "-" + idSampleTumor + "_vs_" + idSampleNormal}

  publishDir directoryMap.snpeff, mode: 'copy'

  input:
    set variantCaller, idPatient, gender, idSampleNormal, idSampleTumor, file(vcf) from vcfForSnpeff
    val snpeffDb from Channel.value(params.genomes[params.genome].snpeffDb)

  output:
    set file("${vcf.baseName}.ann.vcf"), file("${vcf.baseName}_snpEff_genes.txt"), file("${vcf.baseName}_snpEff.csv"), file("${vcf.baseName}_snpEff_summary.html") into snpeffReport

  when: 'snpeff' in tools

  script:
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar \$SNPEFF_HOME/snpEff.jar \
  $snpeffDb \
  -csvStats ${vcf.baseName}_snpEff.csv \
  -v -cancer \
  ${vcf} \
  > ${vcf.baseName}.ann.vcf

  mv snpEff_summary.html ${vcf.baseName}_snpEff_summary.html
  mv ${vcf.baseName}_snpEff.genes.txt ${vcf.baseName}_snpEff_genes.txt
  """
}

verbose ? snpeffReport = snpeffReport.view {"snpEff Reports: $it"} : ''

process RunVEP {
  tag {variantCaller + "-" + idSampleTumor + "_vs_" + idSampleNormal}

  publishDir directoryMap.vep, mode: 'copy'

  input:
    set variantCaller, idPatient, gender, idSampleNormal, idSampleTumor, file(vcf) from vcfForVep

  output:
    set file("${vcf.baseName}_VEP.txt"), file("${vcf.baseName}_VEP.txt_summary.html") into vepReport

  when: 'vep' in tools && variantCaller != 'strelka'

  script:
  if (workflow.profile == 'testing' || workflow.profile == 'docker')
  """
  vep -i $vcf -o ${vcf.baseName}_VEP.txt -offline
  """
  else if (workflow.profile == 'standard' || workflow.profile == 'interactive' || workflow.profile == 'localhost')
  """
  set -euo pipefail
  variant_effect_predictor.pl --cache --dir_cachevariant_effect_predictor.pl --cache --dir_cache /sw/data/uppnex/vep/87/homo_sapiens/87_${params.genome} -o ${vcf.baseName}_VEP.txt
  """
}

verbose ? vepReport = vepReport.view {"VEP Reports: $it"} : ''


process GenerateMultiQCconfig {
  tag {idPatient}

  publishDir directoryMap.multiQC, mode: 'copy'

  input:

  output:
  file("multiqc_config.yaml") into multiQCconfig

  when: 'multiqc' in tools

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
  echo "- Genome: "${params.genome} >> multiqc_config.yaml
  echo "- Step: "${step} >> multiqc_config.yaml
  echo "- Tools: "${tools.join(", ")} >> multiqc_config.yaml
  echo "top_modules:" >> multiqc_config.yaml
  echo "- 'fastqc'" >> multiqc_config.yaml
  echo "- 'picard'" >> multiqc_config.yaml
  echo "- 'samtools'" >> multiqc_config.yaml
  echo "- 'snpeff'" >> multiqc_config.yaml
  """
}

verbose ? multiQCconfig = multiQCconfig.view {"MultiQC config file: $it"} : ''

reportsForMultiQC = Channel.fromPath( 'Reports/{FastQC,MarkDuplicates,SamToolsStats}/*' )
  .mix(bcfReport,
    fastQCreport,
    markDuplicatesReport,
    multiQCconfig,
    recalibratedBamReport,
    snpeffReport,
    vepReport)
  .flatten()
  .unique()
  .toList()

verbose ? reportsForMultiQC = reportsForMultiQC.view {"Reports for MultiQC: $it"} : ''

process RunMultiQC {
  tag {idPatient}

  publishDir directoryMap.multiQC, mode: 'copy'

  input:
    file ('*') from reportsForMultiQC

  output:
    set file("*multiqc_report.html"), file("*multiqc_data") into multiQCReport

    when: 'multiqc' in tools

  script:
  """
  multiqc -f -v .
  """
}

verbose ? multiQCReport = multiQCReport.view {"MultiQC Report: $it"} : ''

/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def checkFile(it) {
  // Check file existence
  final f = file(it)
  if (!f.exists()) {
    exit 1, "Missing file in TSV file: $it, see --help for more information"
  }
  return f
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
  // Loop through all the possible parameters to check their existence and spelling
  return list.every{ checkParameterExistence(it, realList) }
}

def checkParams(it) {
  // Check if params is in this given list
  return it in [
    'callName',
    'call-name',
    'contactMail',
    'contact-mail',
    'genome',
    'genomes',
    'help',
    'project',
    'runTime',
    'run-time',
    'sample',
    'sampleDir',
    'sample-dir',
    'singleCPUMem',
    'single-CPUMem',
    'step',
    'test',
    'tools',
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

def checkStatus(it) {
  // Check if status is correct
  // Status should be only 0 or 1
  // 0 being normal
  // 1 being tumor (or relapse or anything that is not normal...)
  if (!(it in [0, 1])) {
    exit 1, "Status is not recognized in TSV file: $it, see --help for more information"
  }
  return it
}

def checkTSV(it, number) {
  // Check if TSV has the correct number of items in row
  if (it.size() != number) {
    exit 1, "Malformed row in TSV file: $it, see --help for more information"
  }
  return it
}

def checkUppmaxProject() {
  return !((workflow.profile == 'standard' || workflow.profile == 'interactive') && !params.project)
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
    'bcftoolsStats'    : 'Reports/BCFToolsStats',
    'fastQC'           : 'Reports/FastQC',
    'markDuplicatesQC' : 'Reports/MarkDuplicates',
    'multiQC'          : 'Reports/MultiQC',
    'samtoolsStats'    : 'Reports/SamToolsStats',
    'ascat'            : 'VariantCalling/Ascat',
    'freebayes'        : 'VariantCalling/FreeBayes',
    'haplotypecaller'  : 'VariantCalling/HaplotypeCaller',
    'manta'            : 'VariantCalling/Manta',
    'mutect1'          : 'VariantCalling/MuTect1',
    'mutect2'          : 'VariantCalling/MuTect2',
    'strelka'          : 'VariantCalling/Strelka',
    'vardict'          : 'VariantCalling/VarDict',
    'snpeff'           : 'Annotation/SnpEff',
    'vep'              : 'Annotation/VEP'
  ]
}

def defineReferenceMap() {
  if (!(params.genome in params.genomes)) {
    exit 1, "Genome $params.genome not found in configuration"
  }
  genome = params.genomes[params.genome]

  return [
    'acLoci'      : file(genome.acLoci),  // loci file for ascat
    'dbsnp'       : file(genome.dbsnp),
    'dbsnpIndex'  : file(genome.dbsnpIndex),
    'cosmic'      : file(genome.cosmic),  // cosmic VCF with VCF4.1 header
    'cosmicIndex' : file(genome.cosmicIndex),
    'genomeDict'  : file(genome.genomeDict),  // genome reference dictionary
    'genomeFile'  : file(genome.genome),  // FASTA genome reference
    'genomeIndex' : file(genome.genomeIndex),  // genome .fai file
    'bwaIndex'    : file(genome.bwaIndex),  // BWA index files
    // intervals file for spread-and-gather processes (usually chromosome chunks at centromeres)
    'intervals'   : file(genome.intervals),
    // VCFs with known indels (such as 1000 Genomes, Mill’s gold standard)
    'knownIndels' : genome.knownIndels.collect{file(it)},
    'knownIndelsIndex': genome.knownIndelsIndex.collect{file(it)},
  ]
}

def defineStepList() {
  return [
    'preprocessing',
    'realign',
    'recalibrate',
    'skippreprocessing'
  ]
}

def defineToolList() {
  return [
    'ascat',
    'freebayes',
    'haplotypecaller',
    'manta',
    'multiqc',
    'mutect1',
    'mutect2',
    'snpeff',
    'strelka',
    'vardict',
    'vep'
  ]
}

def extractBams(tsvFile) {
  // Channeling the TSV file containing BAM.
  // Format is: "subject gender status sample bam bai"
  return bamFiles = Channel
    .from(tsvFile.readLines())
    .map{line ->
      list      = checkTSV(line.split(),6)
      idPatient = list[0]
      gender    = list[1]
      status    = checkStatus(list[2].toInteger())
      idSample  = list[3]
      bamFile   = checkFile(list[4])
      baiFile   = checkFile(list[5])

      checkFileExtension(bamFile,".bam")
      checkFileExtension(baiFile,".bai")

      [ idPatient, gender, status, idSample, bamFile, baiFile ]
    }
}

def extractFastq(tsvFile) {
  // Channeling the TSV file containing FASTQ.
  // Format is: "subject gender status sample lane fastq1 fastq2"
  return fastqFiles = Channel
    .from(tsvFile.readLines())
    .map{line ->
      list       = checkTSV(line.split(),7)
      idPatient  = list[0]
      gender     = list[1]
      status     = checkStatus(list[2].toInteger())
      idSample   = list[3]
      idRun      = list[4]

      // When testing workflow from github, paths to FASTQ files start from workflow.projectDir and not workflow.launchDir
      fastqFile1 = workflow.commitId && params.test ? checkFile("$workflow.projectDir/${list[5]}") : checkFile("${list[5]}")
      fastqFile2 = workflow.commitId && params.test ? checkFile("$workflow.projectDir/${list[6]}") : checkFile("${list[6]}")

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

  fastq = Channel.create()

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
  return bamFiles = Channel
    .from(tsvFile.readLines())
    .map{line ->
      list       = checkTSV(line.split(),7)
      idPatient  = list[0]
      gender     = list[1]
      status     = checkStatus(list[2].toInteger())
      idSample   = list[3]
      bamFile    = checkFile(list[4])
      baiFile    = checkFile(list[5])
      recalTable = checkFile(list[6])

      checkFileExtension(bamFile,".bam")
      checkFileExtension(baiFile,".bai")
      checkFileExtension(recalTable,".recal.table")

      [ idPatient, gender, status, idSample, bamFile, baiFile, recalTable ]
    }
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
  line = buffered.readLine()
  assert line.startsWith('@')
  line = line.substring(1)
  fields = line.split(' ')[0].split(':')
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

def generateIntervalsForVC(bams, gI) {
  final bamsForVC = Channel.create()
  final vcIntervals = Channel.create()
  (bams, bamsForVC) = bams.into(2)
  (gI, vcIntervals) = gI.into(2)
  bamsForVC = bamsForVC.spread(vcIntervals)
  return [bamsForVC, bams, gI]
}

def grabRevision() {
  return workflow.revision ?: workflow.scriptId.substring(0,10)
}

def helpMessage(version, revision) { // Display help message
  log.info "CANCER ANALYSIS WORKFLOW ~ $version - revision: $revision"
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
  log.info "       Option to configure preprocessing."
  log.info "       Possible values are:"
  log.info "         preprocessing (default, will start workflow with FASTQ files)"
  log.info "         realign (will start workflow with non-realigned BAM files)"
  log.info "         recalibrate (will start workflow with non-recalibrated BAM files)"
  log.info "         skippreprocessing (will start workflow with recalibrated BAM files)"
  log.info "    --tools"
  log.info "       Option to configure which tools to use in the workflow."
  log.info "         Different tools to be separated by commas."
  log.info "       Possible values are:"
  log.info "         mutect1 (use MuTect1 for VC)"
  log.info "         mutect2 (use MuTect2 for VC)"
  log.info "         freebayes (use FreeBayes for VC)"
  log.info "         vardict (use VarDict for VC)"
  log.info "         strelka (use Strelka for VC)"
  log.info "         haplotypecaller (use HaplotypeCaller for normal bams VC)"
  log.info "         manta (use Manta for SV)"
  log.info "         ascat (use Ascat for CNV)"
  log.info "         snpeff (use snpEff for Annotation of Variants)"
  log.info "         vep (use VEP for Annotation of Variants)"
  log.info "    --genome <Genome>"
  log.info "       Use a specific genome version."
  log.info "       Possible values are:"
  log.info "         GRCh37 (Default)"
  log.info "         GRCh38"
  log.info "         smallGRCh37 (Use a small reference (Tests only))"
  log.info "    --help"
  log.info "       you're reading it"
  log.info "    --verbose"
  log.info "       Adds more verbosity to workflow"
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

def retrieveStatus(bamChannel) {
  return bamChannel = bamChannel.map {
    idPatient, gender, bam, bai ->
    tag = bam.baseName.tokenize('.')[0]
    status   = tag[-1..-1].toInteger()
    idSample = tag.take(tag.length()-2)
    [idPatient, gender, status, idSample, bam, bai]
  }
}

def startMessage(version, revision) { // Display start message
  log.info "CANCER ANALYSIS WORKFLOW ~ $version - revision: $revision"
  log.info "Command Line: $workflow.commandLine"
  log.info "Project Dir : $workflow.projectDir"
  log.info "Launch Dir  : $workflow.launchDir"
  log.info "Work Dir    : $workflow.workDir"
  log.info "Genome      : " + params.genome
  log.info "Step        : " + step
  if (tools) {log.info "Tools       : " + tools.join(', ')}
}

def versionMessage(version, revision) { // Display version message
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
  log.info "Genome      : " + params.genome
  log.info "Step        : " + step
  if (tools) {log.info "Tools       : " + tools.join(', ')}
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
