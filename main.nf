#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
=                   C A N C E R    A N A L Y S I S    W O R K F L O W                  =
========================================================================================
 New Cancer Analysis Workflow. Started March 2016.

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

----------------------------------------------------------------------------------------
 Basic command:
 $ nextflow run SciLifeLab/CAW --sample <sample.tsv>

 All variables are configured in the config and sample files. All variables in the config
 file can be reconfigured on the commande line, like:
 --option <option>
----------------------------------------------------------------------------------------
 Workflow processes overview:
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
 - RunVardict - Run VarDict for Variant Calling (Parrallelized processes)
 - ConcatVCF - Merge results from HaplotypeCaller, MuTect1, MuTect2 and VarDict parrallelized processes
 - RunStrelka - Run Strelka for Variant Calling
 - RunManta - Run Manta for Structural Variant Calling
 - RunAlleleCount - Run AlleleCount to prepare for Ascat
 - RunConvertAlleleCounts - Run convertAlleleCounts to prepare for Ascat
 - RunAscat - Run Ascat for CNV
----------------------------------------------------------------------------------------

========================================================================================
=                               C O N F I G U R A T I O N                              =
========================================================================================
*/

revision = grab_git_revision() ?: ''
version  = "v0.9"
refsDefined = true
stepCorrect = true
verbose = false
workflowSteps = []

if( workflow.profile == 'standard' && !params.project ) exit 1, "No UPPMAX project ID found! Use --project" 

/*
 * Arguments handling: Use steps choose which processes to run or to skip
 * Borrowed ideas from https://github.com/guigolab/grape-nf
 * and from https://github.com/NBISweden/wgs-structvar
 */

switch (params) {
  case {params.help} :
    help_message("$version", "$revision")
    exit 1

  case {params.version} :
    version_message("$version", "$revision")
    exit 1

  case {params.steps} : //Getting list of steps from comma-separated strings
    workflowSteps = params.steps.split(',').collect { it.trim() }
    break
}

if (params.verbose) { verbose = true}

/*
 * Define lists of refs, steps and outDir
 */

refs = [
  "genomeFile":  params.genome,       // genome reference
  "genomeIndex": params.genomeIndex,  // genome reference index
  "genomeDict":  params.genomeDict,   // genome reference dictionary
  "kgIndels":    params.kgIndels,     // 1000 Genomes SNPs
  "kgIndex":     params.kgIndex,      // 1000 Genomes SNPs index
  "dbsnp":       params.dbsnp,        // dbSNP
  "dbsnpIndex":  params.dbsnpIndex,   // dbSNP index
  "millsIndels": params.millsIndels,  // Mill's Golden set of SNPs
  "millsIndex":  params.millsIndex,   // Mill's Golden set index
  "cosmic41":    params.cosmic41,     // cosmic vcf file with VCF4.1 header
  "cosmic":      params.cosmic,       // cosmic vcf file
  "intervals":   params.intervals,    // intervals file for spread-and-gather processes (usually chromosome chunks at centromeres)
  "MantaRef":    params.mantaRef,     // copy of the genome reference file
  "MantaIndex":  params.mantaIndex,   // reference index indexed with samtools/0.1.19
  "acLoci":      params.acLoci        // loci file for ascat
]

stepsList = [
  "preprocessing",
  "realign",
  "skipPreprocessing",
  "MuTect1",
  "MuTect2",
  "VarDict",
  "Strelka",
  "HaplotypeCaller",
  "Manta",
  "Ascat"
]

outDir = [
  "preprocessing"   : 'Preprocessing',
  "nonRealigned"    : 'Preprocessing/NonRealigned',
  "recalibrated"    : 'Preprocessing/Recalibrated',
  "VariantCalling"  : 'VariantCalling',
  "MuTect1"         : 'VariantCalling/MuTect1',
  "MuTect2"         : 'VariantCalling/MuTect2',
  "VarDict"         : 'VariantCalling/VarDictJava',
  "Strelka"         : 'VariantCalling/Strelka',
  "HaplotypeCaller" : 'VariantCalling/HaplotypeCaller',
  "Manta"           : 'VariantCalling/Manta',
  "Ascat"           : 'VariantCalling/Ascat'
]

/*
 * Loop through all the references files (params.* defined in the config file) to verify if they really exist.
 */

refs.each { referenceFile, fileToCheck ->
  test = checkRefExistence(referenceFile, fileToCheck)
  !(test) ? refsDefined = false : ""
}

if (!refsDefined) {
  exit 1, 'Missing Reference file(s), see --help for more information'
}

/*
 * Loop through all the possible steps to verify if they are correctly spelled or if they exist.
 */

workflowSteps.each { 
  test = checkStepExistence(it, stepsList)
  !(test) ? stepCorrect = false : ""
}

if (!stepCorrect) {
  exit 1, 'Unknown step parameter(s), see --help for more information'
}

if (('preprocessing' in workflowSteps && ('realign' in workflowSteps || 'skipPreprocessing' in workflowSteps)) || ('realign' in workflowSteps && 'skipPreprocessing' in workflowSteps)) {
  exit 1, 'Please choose only one step between preprocessing, realign and skipPreprocessing, see --help for more information'
}

if (!('preprocessing' in workflowSteps) && !('realign' in workflowSteps) && !('skipPreprocessing' in workflowSteps) ) {
  workflowSteps.add('preprocessing')
}

if (!params.sample) {
  exit 1, 'Missing TSV file, see --help for more information'
}

/*
 * Extract and verify content of TSV file
 */

fastqFiles = Channel.create()

if ('preprocessing' in workflowSteps) {
  fastqFiles = extractFastqFiles(file(params.sample))
  if (verbose) { fastqFiles = fastqFiles.view { "FASTQ files and IDs to process: $it" } }
} else if ( 'realign' in workflowSteps || 'skipPreprocessing' in workflowSteps) {
  bamFiles = extractBamFiles(file(params.sample))
  if (verbose) { bamFiles = bamFiles.view { "Bam files and IDs to process: $it" } }
  fastqFiles.close()
}

start_message("$version", "$revision")

/*
========================================================================================
=                                   P R O C E S S E S                                  =
========================================================================================
*/

process MapReads {
  tag { idRun }
  time { params.runTime * task.attempt }

  input:
  set idPatient, gender, idSample, idRun, file(fq1), file(fq2) from fastqFiles
  file refs["genomeFile"]

  output:
  set idPatient, gender, idSample, idRun, file("${idRun}.bam") into bams

  when: 'preprocessing' in workflowSteps

  script:
  readGroupString="\"@RG\\tID:${idRun}\\tSM:${idSample}\\tLB:${idSample}\\tPL:illumina\""
  """
  #!/bin/bash

  bwa mem -R ${readGroupString} -B 3 -t ${task.cpus} -M ${refs["genomeFile"]} ${fq1} ${fq2} | \
  samtools view -bS -t ${refs["genomeIndex"]} - | \
  samtools sort - > ${idRun}.bam
  """
}

singleBam = Channel.create()
groupedBam = Channel.create()

if ('preprocessing' in workflowSteps) {
  if (verbose) { bams = bams.view { "BAM files before sorting into group or single: $it" } }

  /*
   * Merge or rename bam
   *
   * Borrowed code from https://github.com/guigolab/chip-nf
   * Now, we decide whether bam is standalone or should be merged by idSample (column 2 from channel bams)
   * http://www.nextflow.io/docs/latest/operator.html?highlight=grouptuple#grouptuple
   */

  bams.groupTuple(by:[2])
    .choice(singleBam, groupedBam) { it[3].size() > 1 ? 1 : 0 }
  singleBam = singleBam.map {
    idPatient, gender, idSample, idRun, bam ->
    [idPatient, gender, idSample, bam]
  }
  if (verbose) {
    groupedBam = groupedBam.view { "Grouped BAMs before merge: $it" }
  }
} else {
  singleBam.close()
  groupedBam.close()
}

process MergeBams {
  tag { idSample }

  cpus 1
  queue 'core'
  memory { params.singleCPUMem * task.attempt }
  time { params.runTime * task.attempt }

  input:
  set idPatient, idSample, idRun, file(bam) from groupedBam

  output:
  set idPatient, idSample, file("${idSample}.bam") into mergedBam

  when: 'preprocessing' in workflowSteps

  script:
  """
  #!/bin/bash

  samtools merge ${idSample}.bam ${bam}
  """
}

bamList = Channel.create()

if ('preprocessing' in workflowSteps) {
  /*
   * merge all bams (merged and singles) to a single channel
   */

  if (verbose) {
    singleBam = singleBam.view { "Single BAM: $it" }
    mergedBam = mergedBam.view { "Merged BAM: $it" }
  }

  bamList = mergedBam.mix(singleBam)
  // /!\ It is assumed that every sample of the patient have the same gender
  bamList = bamList.map { idPatient, gender, idSample, bam -> [idPatient[0], gender[0], idSample, bam].flatten() }

  if (verbose) { bamList = bamList.view { "BAM list for MarkDuplicates: $it" } }
} else {
  bamList.close()
}

if ('preprocessing' in workflowSteps || 'realign' in workflowSteps) {
  println file(outDir["preprocessing"]).mkdir() ? "Folder ${outDir['preprocessing']} created" : "Cannot create folder ${outDir['preprocessing']}"
  println file(outDir["nonRealigned"]).mkdir() ? "Folder ${outDir['nonRealigned']} created" : "Cannot create folder ${outDir['nonRealigned']}"
  println file(outDir["recalibrated"]).mkdir() ? "Folder ${outDir['recalibrated']} created" : "Cannot create folder ${outDir['recalibrated']}"
}

process MarkDuplicates {
  tag { idSample }

  publishDir outDir["nonRealigned"], mode: 'copy'

  cpus 1
  queue 'core'
  memory { params.singleCPUMem * task.attempt }
  time { params.runTime * task.attempt }

  input:
  set idPatient, gender, idSample, file(bam) from bamList

  output:
  set idPatient, gender, idSample, file("${idSample}.md.bam"), file("${idSample}.md.bai") into duplicates

  when: 'preprocessing' in workflowSteps

  script:
  """
  #!/bin/bash

  echo -e ${idPatient}\t${gender}\t\$(echo ${idSample} | cut -d_ -f 3)\t\$(echo ${idSample} | cut -d_ -f 1)\t${outDir["nonRealigned"]}/${idSample}.md.bam\t${outDir["nonRealigned"]}/${idSample}.md.bai >> ${workflow.launchDir}/${outDir["nonRealigned"]}/${idPatient}.tsv
  java -Xmx${task.memory.toGiga()}g -jar ${params.picardHome}/MarkDuplicates.jar \
  INPUT=${bam} \
  METRICS_FILE=${bam}.metrics \
  TMP_DIR=. \
  ASSUME_SORTED=true \
  VALIDATION_STRINGENCY=LENIENT \
  CREATE_INDEX=TRUE \
  OUTPUT=${idSample}.md.bam
  """
}

duplicatesInterval = Channel.create()
duplicatesRealign  = Channel.create()

if ('preprocessing' in workflowSteps) {
  duplicatesGrouped  = Channel.create()
  /*
   * create realign intervals, use both tumor+normal as input
   */

  // group the marked duplicates Bams for intervals and realign by overall subject/patient id (idPatient)
  duplicatesGrouped = duplicates.groupTuple(by:[0,1])
} else if ('realign' in workflowSteps) {
  duplicatesGrouped = Channel.create()
  duplicatesGrouped = bamFiles.groupTuple(by:[0,1])
}

if ('preprocessing' in workflowSteps || 'realign' in workflowSteps) {
  // The duplicatesGrouped channel is duplicated, one copy goes to the CreateIntervals process
  // and the other to the IndelRealigner process
  duplicatesGrouped.into(duplicatesInterval, duplicatesRealign)
  if (verbose) {
    duplicatesInterval = duplicatesInterval.view { "BAMs for RealignerTargetCreator grouped by patient ID: $it" }
    duplicatesRealign  = duplicatesRealign.view  { "BAMs for IndelRealigner grouped by patient ID: $it" }
  }
} else {
  duplicatesInterval.close()
  duplicatesRealign.close()
}

process CreateIntervals {
  // Though VCF indexes are not needed explicitly, we are adding them so they will be linked, and not re-created on the fly.
  tag { idPatient }

  time { params.runTime * task.attempt }

  input:
  set idPatient, gender, idSample, file(bam), file(bai) from duplicatesInterval
  file gf from file(refs["genomeFile"])
  file gi from file(refs["genomeIndex"])
  file gd from file(refs["genomeDict"])
  file ki from file(refs["kgIndels"])
  file kix from file(refs["kgIndex"])
  file mi from file(refs["millsIndels"])
  file mix from file(refs["millsIndex"])

  output:
  file("${idPatient}.intervals") into intervals

  when: 'preprocessing' in workflowSteps || 'realign' in workflowSteps

  script:
  bams = bam.collect{"-I $it"}.join(' ')
  """
  #!/bin/bash

  java -Xmx${task.memory.toGiga()}g -jar ${params.gatkHome}/GenomeAnalysisTK.jar \
  -T RealignerTargetCreator \
  $bams \
  -R $gf \
  -known $ki \
  -known $mi \
  -nt ${task.cpus} \
  -XL hs37d5 \
  -o ${idPatient}.intervals
  """
}

if ('preprocessing' in workflowSteps || 'realign' in workflowSteps) {
  if (verbose) { intervals = intervals.view { "Intervals passed to Realign: $it" } }
}

process RealignBams {
  // use nWayOut to split into T/N pair again
  tag { idPatient }

  time { params.runTime * task.attempt }

  input:
  set idPatient, gender, idSample, file(bam), file(bai) from duplicatesRealign
  file gf from file(refs["genomeFile"])
  file gi from file(refs["genomeIndex"])
  file gd from file(refs["genomeDict"])
  file ki from file(refs["kgIndels"])
  file kix from file(refs["kgIndex"])
  file mi from file(refs["millsIndels"])
  file mix from file(refs["millsIndex"])
  file intervals from intervals

  output:
  val(idPatient) into tempIdPatient
  val(gender) into tempGender
  val(idSample) into tempSamples
  file("*.md.real.bam") into tempBams
  file("*.md.real.bai") into tempBais

  when: 'preprocessing' in workflowSteps || 'realign' in workflowSteps

  script:
  bams = bam.collect{"-I $it"}.join(' ')
  """
  #!/bin/bash

  java -Xmx${task.memory.toGiga()}g -jar ${params.gatkHome}/GenomeAnalysisTK.jar \
  -T IndelRealigner \
  $bams \
  -R $gf \
  -targetIntervals $intervals \
  -known $ki \
  -known $mi \
  -XL hs37d5 \
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
      tempSamples.flatten().toSortedList().flatten().merge(
        tempBams.flatten().toSortedList().flatten(),
        tempBais.flatten().toSortedList().flatten()
      ) { sample, bam, bai -> [sample, bam, bai] }
    )
  )
  if (verbose) { realignedBam = realignedBam.view { "realignedBam to BaseRecalibrator: $it" } }
} else {
  realignedBam.close()
}

process CreateRecalibrationTable {
  tag { idSample }

  time { params.runTime * task.attempt }

  input:
  set idPatient, gender, idSample, file(bam), file(bai) from realignedBam
  file refs["genomeFile"]
  file refs["dbsnp"]
  file refs["kgIndels"]
  file refs["millsIndels"]

  output:
  set idPatient, gender, idSample, file(bam), file(bai), file("${idSample}.recal.table") into recalibrationTable

  when: 'preprocessing' in workflowSteps || 'realign' in workflowSteps

  script:
  """
  #!/bin/bash

  java -Xmx${task.memory.toGiga()}g -Djava.io.tmpdir="/tmp" \
  -jar ${params.gatkHome}/GenomeAnalysisTK.jar \
  -T BaseRecalibrator \
  -R ${refs["genomeFile"]} \
  -I $bam \
  -knownSites ${refs["dbsnp"]} \
  -knownSites ${refs["kgIndels"]} \
  -knownSites ${refs["millsIndels"]} \
  -nct ${task.cpus} \
  -XL hs37d5 \
  -l INFO \
  -o ${idSample}.recal.table
  """
}

if ('preprocessing' in workflowSteps || 'realign' in workflowSteps) {
  if (verbose) { recalibrationTable = recalibrationTable.view { "Base recalibrated table for recalibration: $it" } }
}

process RecalibrateBam {
  tag { idSample }

  publishDir outDir["recalibrated"], mode: 'copy'

  // time { params.runTime * task.attempt }

  input:
  set idPatient, gender, idSample, file(bam), file(bai), recalibrationReport from recalibrationTable
  file refs["genomeFile"]

  output:
  set idPatient, gender, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bai") into recalibratedBams

  when: 'preprocessing' in workflowSteps || 'realign' in workflowSteps

  // TODO: ditto as at the previous BaseRecalibrator step, consider using -nct 4
  script:
  """
  #!/bin/bash

  touch ${workflow.launchDir}/${outDir["recalibrated"]}/${idPatient}.tsv
  echo -e ${idPatient}\t${gender}\t\$(echo $idSample | cut -d_ -f 3)\t\$(echo $idSample | cut -d_ -f 1)\t${outDir["recalibrated"]}/${idSample}.recal.bam\t${outDir["recalibrated"]}/${idSample}.recal.bai >> ${workflow.launchDir}/${outDir["recalibrated"]}/${idPatient}.tsv

  java -Xmx${task.memory.toGiga()}g \
  -jar ${params.gatkHome}/GenomeAnalysisTK.jar \
  -T PrintReads \
  -R ${refs["genomeFile"]} \
  -nct ${task.cpus} \
  -I $bam \
  -XL hs37d5 \
  --BQSR $recalibrationReport \
  -o ${idSample}.recal.bam
  """
}

if ('skipPreprocessing' in workflowSteps) {
  recalibratedBams = bamFiles
}

if (verbose) { recalibratedBams = recalibratedBams.view { "Recalibrated Bam for variant Calling: $it" } }

// Here we have a recalibrated bam set, but we need to separate the bam files based on patient status.
// The sample tsv config file which is formatted like: "subject status sample lane fastq1 fastq2"
// cf fastqFiles channel, I decided just to add __status to the sample name to have less changes to do.
// And so I'm sorting the channel if the sample match __0, then it's a normal sample, otherwise tumor.
// Then spread normal over tumor to get each possibilities
// ie. normal vs tumor1, normal vs tumor2, normal vs tumor3
// then copy this channel into channels for each variant calling
// I guess it will still work even if we have multiple normal samples

bamsTumor = Channel.create()
bamsNormal = Channel.create()
bamsAll = Channel.create()
vcfsToMerge = Channel.create()

bamsFHC = Channel.create()
bamsFMT1 = Channel.create()
bamsFMT2 = Channel.create()
bamsFVD = Channel.create()

hcVCF = Channel.create()
mutect1VCF = Channel.create()
mutect2VCF = Channel.create()
vardictVCF = Channel.create()

bamsForStrelka = Channel.create()
bamsForManta = Channel.create()
bamsForAscat = Channel.create()

// separate recalibrate files by filename suffix (which is status): __0 means normal, __1 means tumor recalibrated BAM

recalibratedBams
  .choice(bamsTumor, bamsNormal) { it[2] =~ /__0$/ ? 1 : 0 }

if (verbose) {
  bamsTumor = bamsTumor.view { "Tumor Bam for variant Calling: $it" }
  bamsNormal = bamsNormal.view { "Normal Bam for variant Calling: $it" }
}

// We know that MuTect2 (and other somatic callers) are notoriously slow. To speed them up we are chopping the reference into
// smaller pieces at centromeres (see repeats/centromeres.list), do variant calling by this intervals, and re-merge the VCFs.
// Since we are on a cluster, this can parallelize the variant call processes, and push down the variant call wall clock time significanlty.

// in fact we need two channels: one for the actual genomic region, and an other for names
// without ":", as nextflow is not happy with them (will report as a failed process).
// For region 1:1-2000 the output file name will be something like 1_1-2000_Sample_name.mutect2.vcf
// from the "1:1-2000" string make ["1:1-2000","1_1-2000"]

// define intervals file by --intervals
intervals = Channel.from(file(params.intervals).readLines())
gI = intervals.map{[it,it.replaceFirst(/\:/,"_")]}

if ('HaplotypeCaller' in workflowSteps) {
  (bamsFHC, bamsNormal, gI) = generateIntervalsForVC(bamsNormal, gI)
  if (verbose) { bamsFHC = bamsFHC.view { "Bams with Intervals for HaplotypeCaller: $it" } }
} else {
  bamsFHC.close()
  hcVCF.close()
}

bamsAll = bamsNormal.spread(bamsTumor)

bamsAll = bamsAll.map { // Since idPatientNormal and idPatientTumor are the same, it's removed from bamsAll Channel (same for genderNormal)
  idPatientNormal, genderNormal, idSampleNormal, bamNormal, baiNormal, idPatientTumor, genderTumor, idSampleTumor, bamTumor, baiTumor ->
  [idPatientNormal, genderNormal, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor]
}

if (verbose) { bamsAll = bamsAll.view { "Mapped Recalibrated Bam for variant Calling: $it" } }

if ('MuTect1' in workflowSteps) {
  (bamsFMT1, bamsAll, gI) = generateIntervalsForVC(bamsAll, gI)
  if (verbose) { bamsFMT1 = bamsFMT1.view { "Bams with Intervals for MuTect1: $it" } }
} else {
  bamsFMT1.close()
  mutect1VCF.close()
}

if ('MuTect2' in workflowSteps) {
  (bamsFMT2, bamsAll, gI) = generateIntervalsForVC(bamsAll, gI)
  if (verbose) { bamsFMT2 = bamsFMT2.view { "Bams with Intervals for MuTect2: $it" } }
} else {
  bamsFMT2.close()
  mutect2VCF.close()
}

if ('VarDict' in workflowSteps) {
  (bamsFVD, bamsAll, gI) = generateIntervalsForVC(bamsAll, gI)
  if (verbose) { bamsFVD = bamsFVD.view { "Bams with Intervals for VarDict: $it" } }
} else {
  bamsFVD.close()
  vardictVCF.close()
}

if ('Strelka' in workflowSteps) {
  (bamsAll, bamsForStrelka) = bamsAll.into(2)
  if (verbose) { bamsForStrelka = bamsForStrelka.view { "Bams for Strelka: $it" } }
} else {
  bamsForStrelka.close()
}

if ('Manta' in workflowSteps) {
  (bamsAll, bamsForManta) = bamsAll.into(2)
  if (verbose) { bamsForManta = bamsForManta.view { "Bams for Manta: $it" } }
} else {
  bamsForManta.close()
}

if ('Ascat' in workflowSteps) {
  (bamsAll, bamsForAscat) = bamsAll.into(2)
  if (verbose) { bamsForAscat = bamsForAscat.view { "Bams for Ascat: $it" } }
} else {
  bamsForAscat.close()
}

process RunHaplotypecaller {
  tag { idSampleNormal + "-" + gen_int }

  cpus 1
  queue 'core'
  memory { params.singleCPUMem * task.attempt }
  time { params.runTime * task.attempt }

  input:
  set idPatient, gender, idSampleNormal, file(bamNormal), file(baiNormal), genInt, gen_int from bamsFHC //Are these values `ped to bamNormal already?
  file refs["genomeFile"]
  file refs["genomeIndex"]

  output:
  set val("HaplotypeCaller"), idPatient, gender, idSampleNormal, val("${gen_int}_${idSampleNormal}"), file("${gen_int}_${idSampleNormal}.vcf") into haplotypecallerOutput

  when: 'HaplotypeCaller' in workflowSteps

  //parellelization information: "Many users have reported issues running HaplotypeCaller with the -nct argument, so we recommend using Queue to parallelize HaplotypeCaller instead of multithreading." However, it can take the -nct argument.
  script:
  """
  #!/bin/bash

  java -Xmx${task.memory.toGiga()}g -jar ${params.gatkHome}/GenomeAnalysisTK.jar \
  -T HaplotypeCaller \
  -R ${refs["genomeFile"]} \
  --dbsnp ${refs["dbsnp"]} \
  -I $bamNormal \
  -L \"${genInt}\" \
  -o ${gen_int}_${idSampleNormal}.vcf
  """
}

if ('HaplotypeCaller' in workflowSteps) {
  if (verbose) { haplotypecallerOutput = haplotypecallerOutput.view { "HaplotypeCaller output: $it" } }
  hcVCF = haplotypecallerOutput.map {
    variantCaller, idPatient, gender, idSampleNormal, tag, vcfFile ->
    [variantCaller, idPatient, gender, idSampleNormal, idSampleNormal, vcfFile]
  }.groupTuple(by:[0,1,2,3,4])
  if (verbose) { hcVCF = hcVCF.view { "hcVCF: $it" } }
}

process RunMutect1 {
  tag { idSampleTumor + "-" + gen_int }

  cpus 1 
  queue 'core'
  memory { params.mutect1Mem * task.attempt }
  time { params.runTime * task.attempt }

  input:
  set idPatient, gender, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), genInt, gen_int from bamsFMT1

  output:
  set val("MuTect1"), idPatient, gender, idSampleNormal, idSampleTumor, val("${gen_int}_${idSampleNormal}_${idSampleTumor}"), file("${gen_int}_${idSampleNormal}_${idSampleTumor}.vcf") into mutect1Output

  when: 'MuTect1' in workflowSteps

  script:
  """
  #!/bin/bash

  java -Xmx${task.memory.toGiga()}g -jar ${params.mutect1Home}/muTect.jar \
  -T MuTect \
  -R ${refs["genomeFile"]} \
  --cosmic ${refs["cosmic41"]} \
  --dbsnp ${refs["dbsnp"]} \
  -I:normal $bamNormal \
  -I:tumor $bamTumor \
  -L \"${genInt}\" \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --out ${gen_int}_${idSampleNormal}_${idSampleTumor}.call_stats.out \
  --vcf ${gen_int}_${idSampleNormal}_${idSampleTumor}.vcf
  """
}

if ('MuTect1' in workflowSteps) {
  if (verbose) { mutect1Output = mutect1Output.view { "MuTect1 output: $it" } }
  mutect1VCF = mutect1Output.map {
    variantCaller, idPatient, gender, idSampleNormal, idSampleTumor, tag, vcfFile ->
    [variantCaller, idPatient, gender, idSampleNormal, idSampleTumor, vcfFile]
  }.groupTuple(by:[0,1,2,3,4])
  if (verbose) { mutect1VCF = mutect1VCF.view { "mutect1VCF: $it" } }
}

process RunMutect2 {
  tag { idSampleTumor + "-" + gen_int }

  cpus 1
  queue 'core'
  memory { params.singleCPUMem * task.attempt }
  time { params.runTime * task.attempt }

  input:
  set idPatient, gender, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), genInt, gen_int from bamsFMT2

  output:
  set val("MuTect2"), idPatient, gender, idSampleNormal, idSampleTumor, val("${gen_int}_${idSampleNormal}_${idSampleTumor}"), file("${gen_int}_${idSampleNormal}_${idSampleTumor}.vcf") into mutect2Output

  // we are using MuTect2 shipped in GATK v3.6
  // TODO: the  "-U ALLOW_SEQ_DICT_INCOMPATIBILITY " flag is actually masking a bug in older Picard versions. Using the latest Picard tool
  // this bug should go away and we should _not_ use this flag
  // removed: -nct ${task.cpus} \

  when: 'MuTect2' in workflowSteps

  script:
  """
  #!/bin/bash

  java -Xmx${task.memory.toGiga()}g -jar ${params.gatkHome}/GenomeAnalysisTK.jar \
  -T MuTect2 \
  -R ${refs["genomeFile"]} \
  --cosmic ${refs["cosmic41"]} \
  --dbsnp ${refs["dbsnp"]} \
  -I:normal $bamNormal \
  -I:tumor $bamTumor \
  -U ALLOW_SEQ_DICT_INCOMPATIBILITY \
  -L \"${genInt}\" \
  -o ${gen_int}_${idSampleNormal}_${idSampleTumor}.vcf
  """
}

if ('MuTect2' in workflowSteps) {
  if (verbose) { mutect2Output = mutect2Output.view { "MuTect2 output: $it" } }
  mutect2VCF = mutect2Output.map {
    variantCaller, idPatient, gender, idSampleNormal, idSampleTumor, tag, vcfFile ->
    [variantCaller, idPatient, gender, idSampleNormal, idSampleTumor, vcfFile]
  }.groupTuple(by:[0,1,2,3,4])
  if (verbose) { mutect2VCF = mutect2VCF.view { "mutect2VCF: $it" } }
}

process RunVardict {
  // we are doing the same trick for VarDictJava: running for the whole reference is a PITA, so we are chopping at repeats
  // (or centromeres) where no useful variant calls are expected
  // ~/dev/VarDictJava/build/install/VarDict/bin/VarDict -G /sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/human_g1k_v37_decoy.fasta -f 0.1 -N "tiny" -b "tiny.tumor__1.recal.bam|tiny.normal__0.recal.bam" -z 1 -F 0x500 -c 1 -S 2 -E 3 -g 4 -R "1:131941-141339"
  // we need further filters, but some of the outputs are empty files, confusing the VCF generator script
  tag { idSampleTumor + "-" + gen_int }

  cpus 1
  queue 'core'
  memory { params.singleCPUMem * task.attempt }
  time { params.runTime * task.attempt }

  input:
  set idPatient, gender, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), genInt, gen_int from bamsFVD

  output:
  set val("VarDict"), idPatient, gender, idSampleNormal, idSampleTumor, val("${gen_int}_${idSampleNormal}_${idSampleTumor}"), file("${gen_int}_${idSampleNormal}_${idSampleTumor}.out") into vardictOutput

  when: 'VarDict' in workflowSteps

  script:
  """
  #!/bin/bash

  ${params.vardictHome}/vardict.pl -G ${refs["genomeFile"]} \
  -f 0.01 -N $bamTumor \
  -b "$bamTumor|$bamNormal" \
  -z 1 -F 0x500 \
  -c 1 -S 2 -E 3 -g 4 \
  -R ${genInt} > ${gen_int}_${idSampleNormal}_${idSampleTumor}.out
  """
}

if ('VarDict' in workflowSteps) {
  if (verbose) { vardictOutput = vardictOutput.view { "VarDict output: $it" } }
  vardictVCF = vardictOutput.map {
    variantCaller, idPatient, gender, idSampleNormal, idSampleTumor, tag, vcFile ->
    [variantCaller, idPatient, gender, idSampleNormal, idSampleTumor, vcFile]
  }.groupTuple(by:[0,1,2,3,4])
  if (verbose) { vardictVCF = vardictVCF.view { "vardictVCF: $it" } }
}

if ('HaplotypeCaller' in workflowSteps || 'MuTect1' in workflowSteps || 'MuTect2' in workflowSteps || 'VarDict' in workflowSteps) {
  vcfsToMerge = hcVCF.mix(mutect1VCF, mutect2VCF, vardictVCF)
  if (verbose) { vcfsToMerge = vcfsToMerge.view { "VCFs To be merged: $it" } }
} else {
  vcfsToMerge.close()
}

process ConcatVCF {
  tag { variantCaller == 'HaplotypeCaller' ? variantCaller + "-" + idSampleNormal : variantCaller + "-" + idSampleNormal + "-" + idSampleTumor }

  publishDir "${outDir["VariantCalling"]}/$variantCaller", mode: 'copy'

  input:
  set variantCaller, idPatient, gender, idSampleNormal, idSampleTumor, vcFiles from vcfsToMerge

  output:
  set variantCaller, idPatient, gender, idSampleNormal, idSampleTumor, file("*.vcf") into vcfConcatenated

  when: 'HaplotypeCaller' in workflowSteps || 'MuTect1' in workflowSteps || 'MuTect2' in workflowSteps || 'VarDict' in workflowSteps

  script:
  vcfFiles = vcFiles.join(' -V ')
  if (variantCaller == 'HaplotypeCaller')
    outputFile = "${variantCaller}_${idPatient}_${idSampleNormal}.vcf"
  else 
    outputFile = "${variantCaller}_${idPatient}_${idSampleNormal}_${idSampleTumor}.vcf"

  if (variantCaller == 'VarDict')
    """
    #!/bin/bash

    for i in $vcFiles ;do
      temp=\$(echo \$i | tr -d '[],')
      cat \$temp | ${params.vardictHome}/VarDict/testsomatic.R >> testsomatic.out
    done

    ${params.vardictHome}/VarDict/var2vcf_somatic.pl -f 0.01 -N "${idPatient}_${idSampleNormal}_${idSampleTumor}" testsomatic.out > $outputFile
    """

  else
    """
    #!/bin/bash

    java -Xmx${task.memory.toGiga()}g -cp ${params.gatkHome}/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants --reference ${refs["genomeFile"]} -V $vcfFiles --outputFile $outputFile
    """
}

if ('HaplotypeCaller' in workflowSteps || 'MuTect1' in workflowSteps || 'MuTect2' in workflowSteps || 'VarDict' in workflowSteps) {
  if (verbose) {  vcfConcatenated = vcfConcatenated.view { "vcfConcatenated: $it" } }
}

process RunStrelka {
  tag { idSampleTumor }

  publishDir outDir["Strelka"]

  time { params.runTime * task.attempt }

  input:
  set idPatient, gender, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from bamsForStrelka
  file refs["genomeFile"]
  file refs["genomeIndex"]

  output:
  set val("Strelka"), idPatient, gender, idSampleNormal, idSampleTumor, file("strelka/results/*.vcf") into strelkaOutput

  when: 'Strelka' in workflowSteps

  script:
  """
  #!/bin/bash

  tumorPath=`readlink ${bamTumor}`
  normalPath=`readlink ${bamNormal}`
  ${params.strelkaHome}/bin/configureStrelkaWorkflow.pl \
  --tumor \$tumorPath \
  --normal \$normalPath \
  --ref ${refs["MantaRef"]} \
  --config ${params.strelkaCFG} \
  --output-dir strelka

  cd strelka

  make -j ${task.cpus}
  """
}

if ('Strelka' in workflowSteps) {
  if (verbose) { strelkaOutput = strelkaOutput.view { "Strelka output: $it" } }
}

process RunManta {
  tag { idSampleTumor }

  publishDir outDir["Manta"]

  input:
  set idPatient, gender, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from bamsForManta
  file refs["MantaRef"]
  file refs["MantaIndex"]

  output:
  set val("Manta"), idPatient, gender, idSampleNormal, idSampleTumor, file("${idSampleNormal}_${idSampleTumor}.somaticSV.vcf"),file("${idSampleNormal}_${idSampleTumor}.candidateSV.vcf"),file("${idSampleNormal}_${idSampleTumor}.diploidSV.vcf"),file("${idSampleNormal}_${idSampleTumor}.candidateSmallIndels.vcf") into mantaOutput

  when: 'Manta' in workflowSteps

  //NOTE: Manta is very picky about naming and reference indexes, the input bam should not contain too many _ and the reference index must be generated using a supported samtools version.
  //Moreover, the bam index must be named .bam.bai, otherwise it will not be recognized
  script:
  """
  #!/bin/bash

  mv ${bamNormal} Normal.bam
  mv ${bamTumor} Tumor.bam

  mv ${baiNormal} Normal.bam.bai
  mv ${baiTumor} Tumor.bam.bai

  configManta.py --normalBam Normal.bam --tumorBam Tumor.bam --reference ${refs["MantaRef"]} --runDir MantaDir
  python MantaDir/runWorkflow.py -m local -j ${task.cpus}
  gunzip -c MantaDir/results/variants/somaticSV.vcf.gz > ${idSampleNormal}_${idSampleTumor}.somaticSV.vcf
  gunzip -c MantaDir/results/variants/candidateSV.vcf.gz > ${idSampleNormal}_${idSampleTumor}.candidateSV.vcf
  gunzip -c MantaDir/results/variants/diploidSV.vcf.gz > ${idSampleNormal}_${idSampleTumor}.diploidSV.vcf
  gunzip -c MantaDir/results/variants/candidateSmallIndels.vcf.gz > ${idSampleNormal}_${idSampleTumor}.candidateSmallIndels.vcf
  """
}

if ('Manta' in workflowSteps) {
  if (verbose) { mantaOutput = mantaOutput.view { "Manta output: $it" } }
}

process RunAlleleCount {
  tag {idSampleTumor}
  // Run commands and code from Malin Larsson
  // Based on Jesper Eisfeldt's code

  cpus 1
  queue 'core'
  memory { params.singleCPUMem * task.attempt }

  input:
  set idPatient, gender, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from bamsForAscat
  file refs["genomeFile"]
  file refs["genomeIndex"]
  file refs["acLoci"]

  output:
  set idPatient, gender, idSampleNormal, idSampleTumor, file("${idSampleNormal}.alleleCount"), file("${idSampleTumor}.alleleCount") into alleleCountOutput

  when: 'Ascat' in workflowSteps

  script:
  """
  #!/bin/bash

  alleleCounter -l ${refs["acLoci"]} -r ${refs["genomeFile"]} -b ${bamNormal} -o ${idSampleNormal}.alleleCount;
  alleleCounter -l ${refs["acLoci"]} -r ${refs["genomeFile"]} -b ${bamTumor} -o ${idSampleTumor}.alleleCount;
  """
}

process RunConvertAlleleCounts {
  // R script from Malin Larssons bitbucket repo:
  // https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
  tag {idSampleTumor}

  cpus 1
  queue 'core'
  memory { params.singleCPUMem * task.attempt }

  input:
  set idPatient, gender, idSampleNormal, idSampleTumor, file(normalAlleleCt), file(tumorAlleleCt) from alleleCountOutput

  output:
  set idPatient, gender, idSampleNormal, idSampleTumor, file("${idSampleNormal}.BAF"), file("${idSampleNormal}.LogR"), file("${idSampleTumor}.BAF"), file("${idSampleTumor}.LogR") into convertAlleleCountsOutput

  when: 'Ascat' in workflowSteps

  script:
  """
  #!/bin/bash

  convertAlleleCounts.r ${idSampleTumor} ${tumorAlleleCt} ${idSampleNormal} ${normalAlleleCt} ${gender}
  """
}

process RunAscat {
  // R scripts from Malin Larssons bitbucket repo:
  // https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
  tag {idSampleTumor}

  publishDir outDir["Ascat"]

  cpus 1
  queue 'core'
  memory { params.singleCPUMem * task.attempt }

  input:
  set idPatient, gender, idSampleNormal, idSampleTumor, file(normalBAF), file(normalLogR), file(tumorBAF), file(tumorLogR) from convertAlleleCountsOutput

  output:
  set val("Ascat"), idPatient, gender, idSampleNormal, idSampleTumor, file("${idSampleTumor}.tumour.png"), file("${idSampleTumor}.germline.png"), file("${idSampleTumor}.LogR.PCFed.txt"), file("${idSampleTumor}.BAF.PCFed.txt"), file("${idSampleTumor}.ASPCF.png"), file("${idSampleTumor}.sunrise.png") into ascatOutput

  when: 'Ascat' in workflowSteps

  script:
  """
  #!/bin/env Rscript
  #######################################################################################################
  # Description:
  # R-script for converting output from AlleleCount to BAF and LogR values.
  #
  # Input:
  # AlleleCounter output file for tumor and normal samples
  # The first line should contain a header describing the data
  # The following columns and headers should be present:
  # CHR    POS     Count_A Count_C Count_G Count_T Good_depth
  #
  # Output:
  # BAF and LogR tables (tab delimited text files)
  #######################################################################################################
  source("$baseDir/scripts/ascat.R")
  .libPaths( c( "$baseDir/scripts", .libPaths() ) )
  if(!require(RColorBrewer)){
      source("http://bioconductor.org/biocLite.R")
      biocLite("RColorBrewer", suppressUpdates=TRUE, lib="$baseDir/scripts")
      library(RColorBrewer)
  }
   
  options(bitmapType='cairo')
  tumorbaf = "${tumorBAF}"
  tumorlogr = "${tumorLogR}"
  normalbaf = "${normalBAF}"
  normallogr = "${normalLogR}"
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
  """ //To restore syntaxic coloration: "
}

/*
========================================================================================
=                                   F U N C T I O N S                                  =
========================================================================================
*/

def grab_git_revision() {
  // Borrowed from https://github.com/NBISweden/wgs-structvar

  if ( workflow.commitId ) { // it's run directly from github
    return workflow.commitId.substring(0,10)
  }
  // Try to find the revision directly from git
  head_pointer_file = file("${baseDir}/.git/HEAD")
  if ( ! head_pointer_file.exists() ) {
    return ''
  }
  ref = head_pointer_file.newReader().readLine().tokenize()[1]
  ref_file = file("${baseDir}/.git/$ref")
  if ( ! ref_file.exists() ) {
    return ''
  }
  revision = ref_file.newReader().readLine()
  return revision.substring(0,10)
}

/*
 * Verify parameters and files existence:
 */

def checkRefExistence(referenceFile, fileToCheck) {
  try { assert file(fileToCheck).exists() }
  catch (AssertionError ae) {
    log.info  "Missing references: ${referenceFile} ${fileToCheck}"
    return false
  }
  return true
}

def checkStepExistence(step, stepsList) {
  try { assert stepsList.contains(step) }
  catch (AssertionError ae) {
    println("Unknown parameter: ${step}")
    return false
  }
  return true
}

def checkFileExistence(fileToCheck) {
  try { assert file(fileToCheck).exists() }
  catch (AssertionError ae) {
    exit 1, "Missing file in TSV file: ${fileToCheck}, see --help for more information"
  }
}

def extractFastqFiles(tsvFile) {
  /*
   * Channeling the TSV file containing FASTQ
   * The format is: "subject gender status sample lane fastq1 fastq2"
   * I just added __status to the idSample so that the whole pipeline is still working without having to change anything.
   * I know, it is lazy...
   */
  fastqFiles = Channel
    .from(tsvFile.readLines())
    .map{ line ->
      list        = line.split()
      idPatient   = list[0]
      gender      = list[1]
      idSample    = "${list[3]}__${list[2]}"
      idRun       = list[4]
      fastqFile1  = file(list[5])
      fastqFile2  = file(list[6])

      checkFileExistence(fastqFile1)
      checkFileExistence(fastqFile2)

      [ idPatient, gender, idSample, idRun, fastqFile1, fastqFile2 ]
    }
  return fastqFiles
}

def extractBamFiles(tsvFile) {
  /*
   * Channeling the TSV file containing BAM
   * The format is: "subject gender status sample bam bai"
   * Still with the __status added to the idSample
   */
  bamFiles = Channel
    .from(tsvFile.readLines())
    .map{ line ->
      list        = line.split()
      idPatient   = list[0]
      gender      = list[1]
      idSample    = "${list[3]}__${list[2]}"
      bamFile     = file(list[4])
      baiFile     = file(list[5])

      checkFileExistence(bamFile)
      checkFileExistence(baiFile)

      [ idPatient, gender, idSample, bamFile, baiFile ]
    }
  return bamFiles
}

def generateIntervalsForVC(bams, gI) {
  bamsForVC = Channel.create()
  vcIntervals = Channel.create()
  (bams, bamsForVC) = bams.into(2)
  (gI, vcIntervals) = gI.into(2)
  bamsForVC = bamsForVC.spread(vcIntervals)
  return [bamsForVC, bams, gI]
}

/*
 * Various messages:
 */

def help_message(version, revision) {
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
}

def start_message(version, revision) {
  log.info "CANCER ANALYSIS WORKFLOW ~ $version - revision: $revision"
  log.info "Project     : ${workflow.projectDir}"
  log.info "Directory   : ${workflow.launchDir}"
  log.info "workDir     : ${workflow.workDir}"
  log.info "Steps       : " + workflowSteps.join(", ")
  log.info "Command line: ${workflow.commandLine}"
}

def version_message(version, revision) {
  log.info "CANCER ANALYSIS WORKFLOW"
  log.info "  version $version"
  log.info "  revision: $revision"
  log.info "Git info  : repository - $revision [$workflow.commitId]"
  log.info "Project   : ${workflow.projectDir}"
  log.info "Directory : ${workflow.launchDir}"
}

workflow.onComplete {
  log.info "CANCER ANALYSIS WORKFLOW ~ $version - revision: $revision"
  log.info "Project     : ${workflow.projectDir}"
  log.info "workDir     : ${workflow.workDir}"
  log.info "Command line: ${workflow.commandLine}"
  log.info "Steps       : " + workflowSteps.join(", ")
  log.info "Completed at: ${workflow.complete}"
  log.info "Duration    : ${workflow.duration}"
  log.info "Success     : ${workflow.success}"
  log.info "Exit status : ${workflow.exitStatus}"
  log.info "Error report: ${workflow.errorReport ?: '-'}"
}

workflow.onError {
  log.info "Workflow execution stopped with the following message: ${workflow.errorMessage}"
}