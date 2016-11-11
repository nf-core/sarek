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
 $ nextflow run SciLifeLab/CAW -c <file.config> --sample <sample.tsv>

 All variables are configured in the config and sample files. All variables in the config
 file can be reconfigured on the commande line, like:
 --option <option>
----------------------------------------------------------------------------------------
 Workflow processes overview:
 - MapReads - Map reads
 - MergeBams - Merge BAMs if multilane samples
 - RenameSingleBam - Rename BAM if non-multilane sample
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
 - alleleCount - Run Ascat for CNV
 - convertAlleleCounts - Run Ascat for CNV
 - runASCAT - Run Ascat for CNV
----------------------------------------------------------------------------------------

========================================================================================
=                               C O N F I G U R A T I O N                              =
========================================================================================
*/

revision = grab_git_revision() ?: ''
version  = "v0.8.5"
refsDefined = true
stepCorrect = true
workflowSteps = []
verbose = false

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
  "ascat"
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
  "ascat"           : 'VariantCalling/ascat'
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

workflowSteps.each { step ->
  test = checkStepExistence(step, stepsList)
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
  if (verbose) { fastqFiles = fastqFiles.view { it -> "FASTQ files and IDs to process: $it" } }
} else if ( 'realign' in workflowSteps || 'skipPreprocessing' in workflowSteps) {
  bamFiles = extractBamFiles(file(params.sample))
  if (verbose) { bamFiles = bamFiles.view { it -> "Bam files and IDs to process: $it" } }
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
  file refs["genomeFile"]
  set idPatient, idSample, idRun, file(fq1), file(fq2) from fastqFiles

  output:
  set idPatient, idSample, idRun, file("${idRun}.bam") into bams

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
  if (verbose) { bams = bams.view { it -> "BAM files before sorting into group or single: $it" } }

  /*
   * Merge or rename bam
   *
   * Borrowed code from https://github.com/guigolab/chip-nf
   * Now, we decide whether bam is standalone or should be merged by idSample (column 1 from channel bams)
   * http://www.nextflow.io/docs/latest/operator.html?highlight=grouptuple#grouptuple
   */

  bams.groupTuple(by:[1])
    .choice(singleBam, groupedBam) { it[3].size() > 1 ? 1 : 0 }
  if (verbose) {
    singleBam  = singleBam.view  { it -> "Single BAMs before merge: $it" }
    groupedBam = groupedBam.view { it -> "Grouped BAMs before merge: $it" }
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
  idRun = idRun.sort().join(':')
  """
  #!/bin/bash

  echo -e "idPatient:\t"${idPatient}"\nidSample:\t"${idSample}"\nidRun:\t"${idRun}"\nbam:\t"${bam}"\n" > logInfo
  samtools merge ${idSample}.bam ${bam}
  """
}

process RenameSingleBam {
  // Renaming is totally useless, just to keep file names consistent
  tag { idSample }

  cpus 1
  queue 'core'

  input:
  set idPatient, idSample, idRun, file(bam) from singleBam

  output:
  set idPatient, idSample, file("${idSample}.bam") into singleRenamedBam

  when: 'preprocessing' in workflowSteps

  script:
  """
  #!/bin/bash

  mv ${bam} ${idSample}.bam
  """
}

bamList = Channel.create()

if ('preprocessing' in workflowSteps) {
  /*
   * merge all bams (merged and singles) to a single channel
   */

  if (verbose) {
    singleRenamedBam = singleRenamedBam.view { it -> "Single BAM: $it" }
    mergedBam = mergedBam.view { it -> "Merged BAM: $it" }
  }

  bamList = mergedBam.mix(singleRenamedBam)
  bamList = bamList.map { idPatient, idSample, bam -> [idPatient[0], idSample, bam].flatten() }
  if (verbose) { bamList = bamList.view { it -> "BAM list for MarkDuplicates: $it" } }
} else if ('realign' in workflowSteps) {
  bamList = bamFiles.map { idPatient, idSample, bamFile, baiFile -> [idPatient, idSample, bamFile] }
  if (verbose) { bamList = bamList.view { it -> "BAM list for MarkDuplicates: $it" } }
} else {
  bamList.close()
}

process MarkDuplicates {
  // The output channel is duplicated, one copy goes to the CreateIntervals process
  // and the other to the IndelRealigner process
  tag { idSample }

  cpus 1
  queue 'core'
  memory { params.singleCPUMem * task.attempt }
  time { params.runTime * task.attempt }

  input:
  set idPatient, idSample, file(bam) from bamList

  output:
  set idPatient, idSample, file("${idSample}.md.bam"), file("${idSample}.md.bai") into duplicates

  when: 'preprocessing' in workflowSteps || 'realign' in workflowSteps

  script:
  """
  #!/bin/bash

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

if ('preprocessing' in workflowSteps || 'realign' in workflowSteps) {
  duplicatesGrouped  = Channel.create()
  /*
   * create realign intervals, use both tumor+normal as input
   */

  // group the marked duplicates Bams for intervals and realign by overall subject/patient id (idPatient)
  duplicatesGrouped = duplicates.groupTuple()
  duplicatesGrouped.into(duplicatesInterval, duplicatesRealign)
  if (verbose) {
    duplicatesInterval = duplicatesInterval.view { it -> "BAMs for RealignerTargetCreator grouped by patient ID: $it" }
    duplicatesRealign  = duplicatesRealign.view  { it -> "BAMs for IndelRealigner grouped by patient ID: $it" }
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
  set idPatient, idSample, file(bam), file(bai) from duplicatesInterval
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
  input = mdBam.collect{"-I $it"}.join(' ')
  """
  #!/bin/bash

  java -Xmx${task.memory.toGiga()}g -jar ${params.gatkHome}/GenomeAnalysisTK.jar \
  -T RealignerTargetCreator \
  $input \
  -R $gf \
  -known $ki \
  -known $mi \
  -nt ${task.cpus} \
  -XL hs37d5 \
  -o ${idPatient}.intervals
  """
}

if ('preprocessing' in workflowSteps || 'realign' in workflowSteps) {
  println file(outDir["preprocessing"]).mkdir() ? "Folder ${outDir['preprocessing']} created" : "Cannot create folder ${outDir['preprocessing']}"
  println file(outDir["nonRealigned"]).mkdir() ? "Folder ${outDir['nonRealigned']} created" : "Cannot create folder ${outDir['nonRealigned']}"
  println file(outDir["recalibrated"]).mkdir() ? "Folder ${outDir['recalibrated']} created" : "Cannot create folder ${outDir['recalibrated']}"
  if (verbose) { intervals = intervals.view { it -> "Intervals passed to Realign: $it" } }
}

process RealignBams {
  // use nWayOut to split into T/N pair again
  tag { idPatient }

  publishDir outDir["nonRealigned"], mode: 'copy'

  time { params.runTime * task.attempt }

  input:
  set idPatient, idSample, file(mdBam), file(mdBai) from duplicatesRealign
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
  val(idSample) into tempSamples
  file("*.md.real.bam") into tempBams
  file("*.md.real.bai") into tempBais

  when: 'preprocessing' in workflowSteps || 'realign' in workflowSteps

  script:
  input = mdBam.collect{"-I $it"}.join(' ')

  """
  #!/bin/bash

  rm ${workflow.launchDir}/${outDir["nonRealigned"]}/${idPatient}.tsv
  rm ${workflow.launchDir}/${outDir["recalibrated"]}/${idPatient}.tsv
  touch ${workflow.launchDir}/${outDir["nonRealigned"]}/${idPatient}.tsv

  for i in $idSample ;do
      sample=\$(echo \$i | tr -d '[],')
      echo -e ${idPatient}\t\$(echo \$sample | cut -d_ -f 3)\t\$(echo \$sample | cut -d_ -f 1)\t${outDir["nonRealigned"]}/\$sample.md.real.bam\t${outDir["nonRealigned"]}/\$sample.md.real.bai >> ${workflow.launchDir}/${outDir["nonRealigned"]}/${idPatient}.tsv
  done

  echo -e "idPatient:\t"${idPatient}"\nidSample:\t"${idSample}"\nmdBam:\t"${mdBam}"\nmdBai:\t"${mdBai}"\n" > logInfo
  java -Xmx${task.memory.toGiga()}g -jar ${params.gatkHome}/GenomeAnalysisTK.jar \
  -T IndelRealigner \
  $input \
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
  // We're getting from the Realign process 4 channels (patient, samples bams and bais)
  // So I flatten, sort, and reflatten the samples and the files (bam and bai) channels
  // to get them in the same order (the name of the bam and bai files are based on the sample, so if we sort them they all have the same order ;-))
  // And put them back together, and add the ID patient in the realignedBam channel

  tempSamples = tempSamples.flatten().toSortedList().flatten()
  tempBams = tempBams.flatten().toSortedList().flatten()
  tempBais = tempBais.flatten().toSortedList().flatten()
  tempSamples = tempSamples.merge( tempBams, tempBais ) { s, b, i -> [s, b, i] }
  realignedBam = tempIdPatient.spread(tempSamples)
  if (verbose) { realignedBam = realignedBam.view { it -> "realignedBam to BaseRecalibrator: $it" } }
} else {
  realignedBam.close()
}

process CreateRecalibrationTable {
  tag { idSample }

  time { params.runTime * task.attempt }

  input:
  set idPatient, idSample, file(bam), file(bai) from realignedBam
  file refs["genomeFile"]
  file refs["dbsnp"]
  file refs["kgIndels"]
  file refs["millsIndels"]

  output:
  set idPatient, idSample, file(bam), file(bai), file("${idSample}.recal.table") into recalibrationTable

  when: 'preprocessing' in workflowSteps || 'realign' in workflowSteps

  script:
  """
  #!/bin/bash

  java -Xmx${task.memory.toGiga()}g -Djava.io.tmpdir="/tmp" \
  -jar ${params.gatkHome}/GenomeAnalysisTK.jar \
  -T BaseRecalibrator \
  -R ${refs["genomeFile"]} \
  -I $realignedBamFile \
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
  if (verbose) { recalibrationTable = recalibrationTable.view { it -> "Base recalibrated table for recalibration: $it" } }
}

process RecalibrateBam {
  tag { idSample }

  publishDir outDir["recalibrated"], mode: 'copy'

  // time { params.runTime * task.attempt }

  input:
  set idPatient, idSample, file(bam), file(bai), recalibrationReport from recalibrationTable
  file refs["genomeFile"]

  output:
  set idPatient, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bai") into recalibratedBams

  when: 'preprocessing' in workflowSteps || 'realign' in workflowSteps

  // TODO: ditto as at the previous BaseRecalibrator step, consider using -nct 4
  script:
  """
  #!/bin/bash

  touch ${workflow.launchDir}/${outDir["recalibrated"]}/${idPatient}.tsv
  echo -e ${idPatient}\t\$(echo $idSample | cut -d_ -f 3)\t\$(echo $idSample | cut -d_ -f 1)\t${outDir["recalibrated"]}/${idSample}.recal.bam\t${outDir["recalibrated"]}/${idSample}.recal.bai >> ${workflow.launchDir}/${outDir["recalibrated"]}/${idPatient}.tsv

  java -Xmx${task.memory.toGiga()}g \
  -jar ${params.gatkHome}/GenomeAnalysisTK.jar \
  -T PrintReads \
  -R ${refs["genomeFile"]} \
  -nct ${task.cpus} \
  -I $realignedBamFile \
  -XL hs37d5 \
  --BQSR $recalibrationReport \
  -o ${idSample}.recal.bam
  """
}

if ('skipPreprocessing' in workflowSteps) {
  recalibratedBams = bamFiles
}

if (verbose) { recalibratedBams = recalibratedBams.view { it -> "Recalibrated Bam for variant Calling: $it" } }

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
  .choice(bamsTumor, bamsNormal) { it[1] =~ /__0$/ ? 1 : 0 }

if (verbose) {
  bamsTumor = bamsTumor.view { it -> "Tumor Bam for variant Calling: $it" }
  bamsNormal = bamsNormal.view { it -> "Normal Bam for variant Calling: $it" }
}

// We know that MuTect2 (and other somatic callers) are notoriously slow. To speed them up we are chopping the reference into
// smaller pieces at centromeres (see repeats/centromeres.list), do variant calling by this intervals, and re-merge the VCFs.
// Since we are on a cluster, this can parallelize the variant call processes, and push down the variant call wall clock time significanlty.

// in fact we need two channels: one for the actual genomic region, and an other for names
// without ":", as nextflow is not happy with them (will report as a failed process).
// For region 1:1-2000 the output file name will be something like 1_1-2000_Sample_name.mutect2.vcf
// from the "1:1-2000" string make ["1:1-2000","1_1-2000"]

intervals = Channel // define intervals file by --intervals
  .from(file(params.intervals).readLines())
gI = intervals
  .map {a -> [a,a.replaceFirst(/\:/,"_")]}

if ('HaplotypeCaller' in workflowSteps) {
  bamsForHC = Channel.create()
  hcIntervals = Channel.create()
  (bamsNormal, bamsForHC) = bamsNormal.into(2)
  (gI, hcIntervals) = gI.into(2)
  if (verbose) {
    bamsForHC = bamsForHC.view { it -> "Bams for HaplotypeCaller: $it" }
    hcIntervals = hcIntervals.view { it -> "Intervals for HaplotypeCaller: $it" }
  }
  bamsFHC = bamsForHC.spread(hcIntervals)
  if (verbose) { bamsFHC = bamsFHC.view { it -> "Bams with Intervals for HaplotypeCaller: $it" } }
} else {
  bamsFHC.close()
  hcVCF.close()
}

bamsAll = bamsNormal.spread(bamsTumor)

// Since idPatientNormal and idPatientTumor are the same, it's removed from BamsAll Channel

bamsAll = bamsAll.map {
  idPatientNormal, idSampleNormal, bamNormal, baiNormal, idPatientTumor, idSampleTumor, bamTumor, baiTumor ->
  [idPatientNormal, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor]
}

if (verbose) { bamsAll = bamsAll.view { it -> "Mapped Recalibrated Bam for variant Calling: $it" } }

if ('MuTect1' in workflowSteps) {
  bamsForMutect1 = Channel.create()
  mutect1Intervals = Channel.create()
  (bamsAll, bamsForMutect1) = bamsAll.into(2)
  (gI, mutect1Intervals) = gI.into(2)
  if (verbose) {
    bamsForMutect1 = bamsForMutect1.view { it -> "Bams for MuTect1: $it" }
    mutect1Intervals = mutect1Intervals.view { it -> "Intervals for MuTect1: $it" }
  }
  bamsFMT1 = bamsForMutect1.spread(mutect1Intervals)
  if (verbose) { bamsFMT1 = bamsFMT1.view { it -> "Bams with Intervals for MuTect1: $it" } }
} else {
  bamsFMT1.close()
  mutect1VCF.close()
}

if ('MuTect2' in workflowSteps) {
  bamsForMutect2 = Channel.create()
  mutect2Intervals = Channel.create()
  (bamsAll, bamsForMutect2) = bamsAll.into(2)
  (gI, mutect2Intervals) = gI.into(2)
  if (verbose) {
    bamsForMutect2 = bamsForMutect2.view { it -> "Bams for MuTect2: $it" }
    mutect2Intervals = mutect2Intervals.view { it -> "Intervals for MuTect2: $it" }
  }
  bamsFMT2 = bamsForMutect2.spread(mutect2Intervals)
  if (verbose) { bamsFMT2 = bamsFMT2.view { it -> "Bams with Intervals for MuTect2: $it" } }
} else {
  bamsFMT2.close()
  mutect2VCF.close()
}

if ('VarDict' in workflowSteps) {
  bamsForVardict = Channel.create()
  vardictIntervals = Channel.create()
  (bamsAll, bamsForVardict) = bamsAll.into(2)
  (gI, vardictIntervals) = gI.into(2)
  if (verbose) {
    bamsForVardict = bamsForVardict.view { it -> "Bams for MuTect2: $it" }
    vardictIntervals = vardictIntervals.view { it -> "Intervals for MuTect2: $it" }
  }
  bamsFVD = bamsForVardict.spread(vardictIntervals)
  if (verbose) { bamsFVD = bamsFVD.view { it -> "Bams with Intervals for MuTect2: $it" } }
} else {
  bamsFVD.close()
  vardictVCF.close()
}

if ('Strelka' in workflowSteps) {
  (bamsAll, bamsForStrelka) = bamsAll.into(2)
  if (verbose) { bamsForStrelka = bamsForStrelka.view { it -> "Bams for Strelka: $it" } }
} else {
  bamsForStrelka.close()
}

if ('Manta' in workflowSteps) {
  (bamsAll, bamsForManta) = bamsAll.into(2)
  if (verbose) { bamsForManta = bamsForManta.view { it -> "Bams for Manta: $it" } }
} else {
  bamsForManta.close()
}

if ('ascat' in workflowSteps) {
  (bamsAll, bamsForAscat) = bamsAll.into(2)
  if (verbose) { bamsForAscat = bamsForAscat.view { it -> "Bams for ascat: $it" } }
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
  set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), genInt, gen_int from bamsFHC //Are these values `ped to bamNormal already?
  file refs["genomeFile"]
  file refs["genomeIndex"]

  output:
  set idPatient, idSampleNormal, val("${gen_int}_${idSampleNormal}"), file("${gen_int}_${idSampleNormal}.vcf") into haplotypeCallerVariantCallingOutput

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
  if (verbose) { haplotypeCallerVariantCallingOutput = haplotypeCallerVariantCallingOutput.view { it -> "HaplotypeCaller output: $it" } }
  haplotypeCallerVariantCallingOutput = haplotypeCallerVariantCallingOutput.map {
    idPatient, idSampleNormal, tag, vcfFile ->
    [idPatient, idSampleNormal, idSampleNormal, vcfFile]
  }
  hcVCF = Channel.from('HaplotypeCaller').spread(haplotypeCallerVariantCallingOutput)
  hcVCF = hcVCF.groupTuple(by:[0,1,2,3])
  if (verbose) { hcVCF = hcVCF.view { it -> "hcVCF: $it" } }
}

process RunMutect1 {
  tag { idSampleTumor + "-" + gen_int }

  cpus 1 
  queue 'core'
  memory { params.mutect1Mem * task.attempt }
  time { params.runTime * task.attempt }

  input:
  set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), genInt, gen_int from bamsFMT1

  output:
  set idPatient, idSampleNormal, idSampleTumor, val("${gen_int}_${idSampleNormal}_${idSampleTumor}"), file("${gen_int}_${idSampleNormal}_${idSampleTumor}.vcf") into mutect1VariantCallingOutput

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
  if (verbose) { mutect1VariantCallingOutput = mutect1VariantCallingOutput.view { it -> "MuTect1 output: $it" } }
  mutect1VariantCallingOutput = mutect1VariantCallingOutput.map {
    idPatient, idSampleNormal, idSampleTumor, tag, vcfFile ->
    [idPatient, idSampleNormal, idSampleTumor, vcfFile]
  }
  mutect1VCF = Channel.from('MuTect1').spread(mutect1VariantCallingOutput)
  mutect1VCF = mutect1VCF.groupTuple(by:[0,1,2,3])
  if (verbose) { mutect1VCF = mutect1VCF.view { it -> "mutect1VCF: $it" } }
}

process RunMutect2 {
  tag { idSampleTumor + "-" + gen_int }

  cpus 1
  queue 'core'
  memory { params.singleCPUMem * task.attempt }
  time { params.runTime * task.attempt }

  input:
  set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), genInt, gen_int from bamsFMT2

  output:
  set idPatient, idSampleNormal, idSampleTumor, val("${gen_int}_${idSampleNormal}_${idSampleTumor}"), file("${gen_int}_${idSampleNormal}_${idSampleTumor}.vcf") into mutect2VariantCallingOutput

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
  if (verbose) { mutect2VariantCallingOutput = mutect2VariantCallingOutput.view { it -> "MuTect2 output: $it" } }
  mutect2VariantCallingOutput = mutect2VariantCallingOutput.map {
    idPatient, idSampleNormal, idSampleTumor, tag, vcfFile ->
    [idPatient, idSampleNormal, idSampleTumor, vcfFile]
  }
  mutect2VCF = Channel.from('MuTect2').spread(mutect2VariantCallingOutput)
  mutect2VCF = mutect2VCF.groupTuple(by:[0,1,2,3])
  if (verbose) { mutect2VCF = mutect2VCF.view { it -> "mutect2VCF: $it" } }
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
  set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), genInt, gen_int from bamsFVD

  output:
  set idPatient, idSampleNormal, idSampleTumor, val("${gen_int}_${idSampleNormal}_${idSampleTumor}"), file("${gen_int}_${idSampleNormal}_${idSampleTumor}.out") into vardictVariantCallingOutput

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
  if (verbose) { vardictVariantCallingOutput = vardictVariantCallingOutput.view { it -> "VarDict output: $it" } }
  vardictVariantCallingOutput = vardictVariantCallingOutput.map {
    idPatient, idSampleNormal, idSampleTumor, tag, vcfFile ->
    [idPatient, idSampleNormal, idSampleTumor, vcfFile]
  }
  vardictVCF = Channel.from('VarDict').spread(vardictVariantCallingOutput)
  vardictVCF = vardictVCF.groupTuple(by:[0,1,2,3])
  if (verbose) { vardictVCF = vardictVCF.view { it -> "vardictVCF: $it" } }
}

if ('HaplotypeCaller' in workflowSteps || 'MuTect1' in workflowSteps || 'MuTect2' in workflowSteps || 'VarDict' in workflowSteps) {
  vcfsToMerge = hcVCF.mix(mutect1VCF, mutect2VCF, vardictVCF)
  if (verbose) { vcfsToMerge = vcfsToMerge.view { it -> "VCFs To be merged: $it" } }
} else {
  vcfsToMerge.close()
}

process ConcatVCF {
  tag { variantCaller == 'HaplotypeCaller' ? variantCaller + "-" + idSampleNormal : variantCaller + "-" + idSampleNormal + "-" + idSampleTumor }

  publishDir "${outDir["VariantCalling"]}/$variantCaller", mode: 'copy'

  input:
  set variantCaller, idPatient, idSampleNormal, idSampleTumor, vcFiles from vcfsToMerge

  output:
  file "*.vcf" into vcfConcatenated

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
      cat \$temp | ${params.vardictHome}/testsomatic.R >> testsomatic.out
    done

    ${params.vardictHome}/var2vcf_somatic.pl -f 0.01 -N "${idPatient}_${idSampleNormal}_${idSampleTumor}" testsomatic.out > outputFile
    """

  else
    """
    #!/bin/bash

    java -Xmx${task.memory.toGiga()}g -cp ${params.gatkHome}/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants --reference ${refs["genomeFile"]}  -V $vcfFiles --outputFile $outputFile
    """
}

if ('HaplotypeCaller' in workflowSteps || 'MuTect1' in workflowSteps || 'MuTect2' in workflowSteps || 'VarDict' in workflowSteps) {
  if (verbose) {  vcfConcatenated = vcfConcatenated.view { it -> "vcfConcatenated: $it" } }
}

process RunStrelka {
  tag { idSampleTumor }

  publishDir outDir["Strelka"]

  time { params.runTime * task.attempt }

  input:
  set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from bamsForStrelka
  file refs["genomeFile"]
  file refs["genomeIndex"]

  output:
  set idPatient, val("${gen_int}_${idSampleNormal}_${idSampleTumor}"), file("strelka/results/*.vcf") into strelkaVariantCallingOutput

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
  if (verbose) { strelkaVariantCallingOutput = strelkaVariantCallingOutput.view { it -> "Strelka output: $it" } }
}

process Manta {
  tag { idSampleTumor }

  publishDir outDir["Manta"]

  input:
    file refs["MantaRef"]
    file refs["MantaIndex"]
    set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from bamsForManta

  output:
   set idPatient, val("${idSampleNormal}_${idSampleTumor}"),file("${idSampleNormal}_${idSampleTumor}.somaticSV.vcf"),file("${idSampleNormal}_${idSampleTumor}.candidateSV.vcf"),file("${idSampleNormal}_${idSampleTumor}.diploidSV.vcf"),file("${idSampleNormal}_${idSampleTumor}.candidateSmallIndels.vcf")  into mantaVariantCallingOutput

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
  if (verbose) { mantaVariantCallingOutput = mantaVariantCallingOutput.view { it -> "Manta output: $it" } }
}

process alleleCount{
  tag {idSampleTumor}
  //  #!/usr/bin/env nextflow
  // This module runs ascat preprocessing and run
  // Commands and code from Malin Larsson
  // Module based on Jesper Eisfeldt's code

  /* Workflow:
      First: run alleleCount on both normal and tumor (each its own process
      Second: Run R script to process allele counts into logR and BAF values
      Third: run ascat
  */

  // 1)
  // module load bioinfo-tools
  // module load alleleCount
  // alleleCounter -l /sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/1000G_phase3_20130502_SNP_maf0.3.loci -r /sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/human_g1k_v37_decoy.fasta -b sample.bam -o sample.allecount

  cpus 1
  queue 'core'
  memory { params.singleCPUMem * task.attempt }

  input:
  file refs["genomeFile"]
  file refs["genomeIndex"]
  file refs["acLoci"]
  // file normal_bam
  set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from bamsForAscat

  output:
  set idPatient, idSampleNormal, idSampleTumor, file("${idSampleNormal}.alleleCount"), file("${idSampleTumor}.alleleCount") into allele_count_output

  when: 'ascat' in workflowSteps

  script:
  """
  #!/bin/bash

  alleleCounter -l ${refs["acLoci"]} -r ${refs["genomeFile"]} -b ${bamNormal} -o ${idSampleNormal}.alleleCount;
  alleleCounter -l ${refs["acLoci"]} -r ${refs["genomeFile"]} -b ${bamTumor} -o ${idSampleTumor}.alleleCount;
  """
} // end process alleleCount

process convertAlleleCounts {
  tag {idSampleTumor}
  // ascat step 2/3
  // converte allele counts
  // R script from Malin Larssons bitbucket repo:
  // https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
  //
  // copyright?

  // prototype: "Rscript convertAlleleCounts.r tumorid tumorac normalid normalac gender"

  cpus 1
  queue 'core'
  memory { params.singleCPUMem * task.attempt }

  input:
  set idPatient, idSampleNormal, idSampleTumor, file(normalAlleleCt), file(tumorAlleleCt) from allele_count_output
  //file ${refs["scriptDir"]}/convertAlleleCounts.r

  output:
  set idPatient, idSampleNormal, idSampleTumor, file("${idSampleNormal}.BAF"), file("${idSampleNormal}.LogR"), file("${idSampleTumor}.BAF"), file("${idSampleTumor}.LogR") into convert_ac_output

  when: 'ascat' in workflowSteps

  script:
  """
  #!/bin/bash

  convertAlleleCounts.r ${idSampleTumor} ${tumorAlleleCt} ${idSampleNormal} ${normalAlleleCt} ${refs["gender"]}
  """


} // end process convertAlleleCounts

process runASCAT {
  tag {idSampleTumor}
  // ascat step 3/3
  // run ascat
  // R scripts from Malin Larssons bitbucket repo:
  // https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
  //
  // copyright?
  //
  // prototype: "Rscript run_ascat.r tumor_baf tumor_logr normal_baf normal_logr"

  cpus 1
  queue 'core'
  memory { params.singleCPUMem * task.attempt }

  input:
  set idPatient, idSampleNormal, idSampleTumor, file(normalBAF), file(normalLogR), file(tumorBAF), file(tumorLogR) from convert_ac_output

  output:
  file "ascat.done"

  when: 'ascat' in workflowSteps

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
  str(ascat.output)
  plot(sort(ascat.output\$aberrantcellfraction))
  plot(density(ascat.output\$ploidy))
  """
  // the following works when ascat.R is in the run_ascat.r file and the run_ascat.r file is in bin/
  //  run_ascat.r ${tumorBAF} ${tumorLogR} ${normalBAF} ${normalLogR}
  //  touch ascat.done

} // end process runASCAT

/*
 * add process for convert allele counts
 * add process for runASCAT.r
*/


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
   * The format is: "subject status sample lane fastq1 fastq2"
   * I just added __status to the idSample so that the whole pipeline is still working without having to change anything.
   * I know, it is lazy...
   */
  fastqFiles = Channel
    .from(tsvFile.readLines())
    .map{ line ->
      list        = line.split()
      idPatient   = list[0]
      idSample    = "${list[2]}__${list[1]}"
      idRun       = list[3]
      fastqFile1  = file(list[4])
      fastqFile2  = file(list[5])

      checkFileExistence(fastqFile1)
      checkFileExistence(fastqFile2)

      [ idPatient, idSample, idRun, fastqFile1, fastqFile2 ]
    }
  return fastqFiles
}

def extractBamFiles(tsvFile) {
  /*
   * Channeling the TSV file containing BAM
   * The format is: "subject status sample bam bai"
   * Still with the __status added to the idSample
   */
  bamFiles = Channel
    .from(tsvFile.readLines())
    .map{ line ->
      list        = line.split()
      idPatient   = list[0]
      idSample    = "${list[2]}__${list[1]}"
      bamFile     = file(list[3])
      baiFile     = file(list[4])

      checkFileExistence(bamFile)
      checkFileExistence(baiFile)

      [ idPatient, idSample, bamFile, baiFile ]
    }
  return bamFiles
}

/*
 * Various messages:
 */

def help_message(version, revision) {
  log.info "CANCER ANALYSIS WORKFLOW ~ $version - revision: $revision"
  log.info "    Usage:"
  log.info "       nextflow run SciLifeLab/CAW -c <file.config> --sample <sample.tsv> [--steps STEP[,STEP]]"
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
  log.info "         ascat (use ascat for CNV)"
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