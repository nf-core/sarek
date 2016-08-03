#!/usr/bin/env nextflow

/*
========================================================================================
=                   C A N C E R    A N A L Y S I S    W O R K F L O W                  =
========================================================================================
 New Cancer Analysis Workflow. Started March 2016.

 @Authors
 Pelin Akan <pelin.akan@scilifelab.se>
 Jesper Eisfeldt <jesper.eisfeldt@scilifelab.se>
 Maxime Garcia <maxime.garcia@scilifelab.se>
 Szilveszter Juhos <szilveszter.juhos@scilifelab.se>
 Max Käller <max.kaller@scilifelab.se>
 Malin Larsson <malin.larsson@scilifelab.se>
 Björn Nystedt <bjorn.nystedt@scilifelab.se>
 Pall Olason <pall.olason@scilifelab.se>
----------------------------------------------------------------------------------------
@Licence
 The MIT License (MIT)

Copyright (c) 2016 SciLifeLab

Permission is hereby granted, free of charge, to any person obtaining a copy of this
software and associated documentation files (the "Software"), to deal in the Software
without restriction, including without limitation the rights to use, copy, modify, merge,
publish, distribute, sublicense, and/or sell copies of the Software, and to permit
persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
----------------------------------------------------------------------------------------
 Basic command:
 $ nextflow run MultiFQtoVC.nf -c <file.config> --sample <sample.tsv>

 All variables are configured in the config and sample files. All variables in the config
 file can be reconfigured on the commande line, like:
 --option <option>
----------------------------------------------------------------------------------------
 Workflow process overview:
 - Mapping - Map reads with BWA
 - MergeBam - Merge BAMs if multilane samples
 - RenameSingleBam - Rename BAM if non-multilane sample
 - MarkDuplicates - using Picard
 - CreateIntervals - using GATK
 - Realign - using GATK
 - CreateRecalibrationTable - using GATK
 - RecalibrateBam - using GATK
 - RunMutect2 - using MuTect2 shipped in GATK v3.6
 - VarDict - run VarDict on multiple intervals
 - VarDictCollatedVCF - merge Vardict result
----------------------------------------------------------------------------------------

========================================================================================
=                               C O N F I G U R A T I O N                              =
========================================================================================
*/

String version = "0.0.2"
String dateUpdate = "2016-08-01"

/*
 * Get some basic informations about the workflow
 * to get more informations use -with-trace or -with-timeline
 */

workflow.onComplete {
  text = Channel.from(
    "CANCER ANALYSIS WORKFLOW",
    "Version     : $version",
    "Command line: ${workflow.commandLine}",
    "Completed at: ${workflow.complete}",
    "Duration    : ${workflow.duration}",
    "Success     : ${workflow.success}",
    "workDir     : ${workflow.workDir}",
    "Exit status : ${workflow.exitStatus}",
    "Error report: ${workflow.errorReport ?: '-'}")
  text.subscribe { log.info "$it" }
}

/*
 * Basic argument handling
 */

switch (params) {
  case {params.help} :
    text = Channel.from(
      "CANCER ANALYSIS WORKFLOW ~ version $version",
      "    Usage:",
      "       nextflow run MultiFQtoVC.nf -c <file.config> --sample <sample.tsv>",
      "    --help",
      "       you're reading it",
      "    --version",
      "       displays version number")
    text.subscribe { println "$it" }
    exit 1

  case {params.version} :
    text = Channel.from(
      "CANCER ANALYSIS WORKFLOW",
      "  Version $version",
      "  Last update on $dateUpdate",
      "Project : $workflow.projectDir",
      "Cmd line: $workflow.commandLine")
    text.subscribe { println "$it" }
    exit 1
}

/*
 * Use this closure to loop through all the parameters.
 * We can get an AssertionError exception from the file() method as well.
 */

parametersDefined = true
CheckExistence = {
  referenceFile, fileToCheck ->
  try {
    referenceFile = file(fileToCheck)
    assert referenceFile.exists()
  }
  catch (AssertionError ae) {
    println("Missing file: ${referenceFile} ${fileToCheck}")
    parametersDefined = false;
  }
}

refs = [
  "genomeFile":   params.genome,      // genome reference
  "genomeIndex":  params.genomeIndex, // genome reference index
  "genomeDict":   params.genomeDict,  // genome reference dictionary
  "kgIndels":     params.kgIndels,    // 1000 Genomes SNPs
  "kgIndex":      params.kgIndex,     // 1000 Genomes SNPs index
  "dbsnp":        params.dbsnp,       // dbSNP
  "dbsnpIndex":   params.dbsnpIndex,  // dbSNP index
  "millsIndels":  params.millsIndels, // Mill's Golden set of SNPs
  "millsIndex":   params.millsIndex,  // Mill's Golden set index
  "sample":       params.sample,      // the sample sheet (multilane data refrence table, see below)
  "cosmic":       params.cosmic       // cosmic vcf file
]

refs.each(CheckExistence)

if (!parametersDefined) {
  text = Channel.from(
    "CANCER ANALYSIS WORKFLOW ~ version $version",
    "Missing file or parameter: please review your config file.",
    "    Usage",
    "       nextflow run MultiFQtoVC.nf -c <file.config> --sample <sample.tsv>")
  text.subscribe { println "$it" }
  exit 1
}

/*
 * Time to check the sample file. Its format is like: "subject status sample lane fastq1 fastq2":
HCC1954 0 HCC1954.blood HCC1954.blood_1 data/HCC1954.normal_S22_L001_R1_001.fastq.gz data/HCC1954.normal_S22_L001_R2_001.fastq.gz
HCC1954 0 HCC1954.blood HCC1954.blood_2 data/HCC1954.normal_S22_L002_R1_001.fastq.gz data/HCC1954.normal_S22_L002_R2_001.fastq.gz
HCC1954 1 HCC1954.tumor HCC1954.tumor_1 data/HCC1954.tumor1_S23_L001_R1_001.fastq.gz data/HCC1954.tumor1_S23_L001_R2_001.fastq.gz
HCC1954 1 HCC1954.tumor HCC1954.tumor_2 data/HCC1954.tumor1_S23_L002_R1_001.fastq.gz data/HCC1954.tumor1_S23_L002_R2_001.fastq.gz
HCC1954 1 HCC1954.relapse HCC1954.relapse_1 data/HCC1954.tumor2_S24_L001_R1_001.fastq.gz data/HCC1954.tumor2_S24_L001_R2_001.fastq.gz
HCC1954 1 HCC1954.relapse HCC1954.relapse_2 data/HCC1954.tumor2_S24_L002_R1_001.fastq.gz data/HCC1954.tumor2_S24_L002_R2_001.fastq.gz
HCC1954 1 HCC1954.9746123 HCC1954.9746123_1 data/HCC1954.tumor3_S25_L001_R1_001.fastq.gz data/HCC1954.tumor3_S25_L001_R2_001.fastq.gz
HCC1954 1 HCC1954.9746123 HCC1954.9746123_2 data/HCC1954.tumor3_S25_L002_R1_001.fastq.gz data/HCC1954.tumor3_S25_L002_R2_001.fastq.gz
 */

sampleTSVconfig = file(params.sample)

if (!params.sample) {
  text = Channel.from(
    "CANCER ANALYSIS WORKFLOW ~ version $version",
    "Missing the sample TSV config file: please specify it.",
    "    Usage",
    "       nextflow run MultiFQtoVC.nf -c <file.config> --sample <sample.tsv>")
  text.subscribe { println "$it" }
  exit 1
}

/*
 * Read config file, it's "subject status sample lane fastq1 fastq2"
 * let's channel this out for mapping
 */

// [maxime] I just added __status to the sample ID so that the whole pipeline is still working without having to change anything.
// [maxime] I know, it is lazy...

fastqFiles = Channel
  .from(sampleTSVconfig.readLines())
  .map {line ->
    list        = line.split()
    idPatient   = list[0]
    idSample    = "${list[2]}__${list[1]}"
    idRun       = list[3]
    fastqFile1  = file(list[4])
    fastqFile2  = file(list[5])
    [ idPatient, idSample, idRun, fastqFile1, fastqFile2 ]
}

/*
========================================================================================
=                                   P R O C E S S E S                                  =
========================================================================================
*/

fastqFiles = logChannelContent("FASTQ files and IDs to process: ", fastqFiles)

process Mapping {

  module 'bioinfo-tools'
  module 'bwa/0.7.8'
  module 'samtools/1.3'

  memory { 16.GB * task.attempt }
  time { 20.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
  maxRetries 3
  maxErrors '-1'
  cpus 8

  input:
  file refs["genomeFile"]
  set idPatient, idSample, idRun, file(fq1), file(fq2) from fastqFiles

  output:
  set idPatient, idSample, idRun, file("${idRun}.bam") into bams

  // here I use params.genome for bwa ref so I don't have to link to all bwa index files

  script:
  readGroupString="\"@RG\\tID:${idRun}\\tSM:${idSample}\\tLB:${idSample}\\tPL:illumina\""

  """
  bwa mem -R ${readGroupString} -B 3 -t ${task.cpus} -M ${refs["genomeFile"]} ${fq1} ${fq2} | \
  samtools view -bS -t ${refs["genomeIndex"]} - | \
  samtools sort - > ${idRun}.bam
  """
}

bams  = logChannelContent("BAM files before sorting into group or single:", bams)

/*
 * Borrowed code from chip.nf (https://github.com/guigolab/chip-nf)
 *
 * Now, we decide whether bam is standalone or should be merged by sample (id (column 1) from channel bams)
 * http://www.nextflow.io/docs/latest/operator.html?highlight=grouptuple#grouptuple
 *
 */

// Merge or rename bam

singleBam  = Channel.create()
groupedBam = Channel.create()

bams.groupTuple(by:[1])
  .choice(singleBam, groupedBam) { it[3].size() > 1 ? 1 : 0 }

singleBam  = logChannelContent("Single BAMs before merge:", singleBam)
groupedBam = logChannelContent("Grouped BAMs before merge:", groupedBam)

process MergeBam {

  module 'bioinfo-tools'
  module 'samtools/1.3'

  memory { 8.GB * task.attempt }
  time { 16.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
  maxRetries 3
  maxErrors '-1'

  input:
  set idPatient, idSample, idRun, file(bam) from groupedBam

  output:
  set idPatient, idSample, idRun, file("${idSample}.bam") into mergedBam

  script:
  idRun = idRun.sort().join(':')

  """
  echo -e "idPatient:\t"${idPatient}"\nidSample:\t"${idSample}"\nidRun:\t"${idRun}"\nbam:\t"${bam}"\n" > logInfo
  samtools merge ${idSample}.bam ${bam}
  """
}

// [maxime] Renaming is totally useless, but it is more consistent with the rest of the pipeline

process RenameSingleBam {

  input:
  set idPatient, idSample, idRun, file(bam) from singleBam

  output:
  set idPatient, idSample, idRun, file("${idSample}.bam") into singleRenamedBam

  script:
  idRun = idRun.sort().join(':')

  """
  mv ${bam} ${idSample}.bam
  """
}

singleRenamedBam = logChannelContent("SINGLES: ", singleRenamedBam)
mergedBam = logChannelContent("GROUPED: ", mergedBam)

/*
 * merge all bams (merged and singles) to a single channel
 */

bamList = Channel.create()
bamList = mergedBam.mix(singleRenamedBam)
bamList = bamList.map { idPatient, idSample, idRun, bam -> [idPatient[0], idSample, bam].flatten() }

bamList = logChannelContent("BAM list for MarkDuplicates: ",bamList)

/*
 *  mark duplicates all bams
 */

process MarkDuplicates {

  memory { 16.GB * task.attempt }
  time { 16.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
  maxRetries 3
  maxErrors '-1'

  module 'bioinfo-tools'
  module 'picard/1.118'

  input:
  set idPatient, idSample, file(bam) from bamList

  // Channel content should be in the log before
  // The output channels are duplicated nevertheless, one copy goes to RealignerTargetCreator (CreateIntervals)
  // and the other to IndelRealigner

  output:
  set idPatient, idSample, file("${idSample}.md.bam"), file("${idSample}.md.bai") into duplicatesForInterval
  set idPatient, idSample, file("${idSample}.md.bam"), file("${idSample}.md.bai") into duplicatesForRealignement

  """
  echo -e "idPatient:\t"${idPatient}"\nidSample:\t"${idSample}"\nbam:\t"${bam}"\n" > logInfo
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

/*
 * create realign intervals, use both tumor+normal as input
 */

duplicatesForInterval = logChannelContent("BAMs for IndelRealigner before groupTuple: ", duplicatesForInterval)

// group the marked duplicates Bams intervals by overall subject/patient id (idPatient)
duplicatesInterval = Channel.create()
duplicatesInterval = duplicatesForInterval.groupTuple()
duplicatesInterval = logChannelContent("BAMs for RealignerTargetCreator grouped by overall subject/patient ID: ", duplicatesInterval)

duplicatesForRealignement = logChannelContent("BAMs for IndelRealigner before groupTuple: ",  duplicatesForRealignement)

// group the marked duplicates Bams for realign by overall subject/patient id (idPatient)
duplicatesRealign  = Channel.create()
duplicatesRealign  = duplicatesForRealignement.groupTuple()
duplicatesRealign  = logChannelContent("BAMs for IndelRealigner grouped by overall subject/patient ID: ", duplicatesRealign)

/*
 * Creating target intervals for indel realigner.
 * Though VCF indexes are not needed explicitly, we are adding them so they will be linked, and not re-created on the fly.
 */

process CreateIntervals {

  cpus 8
  memory { 8.GB * task.attempt }
  time { 8.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
  maxRetries 3
  maxErrors '-1'

  input:
  set idPatient, idSample, file(mdBam), file(mdBai) from duplicatesInterval
  file gf from file(refs["genomeFile"])
  file gi from file(refs["genomeIndex"])
  file gd from file(refs["genomeDict"])
  file ki from file(refs["kgIndels"])
  file kix from file(refs["kgIndex"])
  file mi from file(refs["millsIndels"])
  file mix from file(refs["millsIndex"])

  output:
  file("${idPatient}.intervals") into intervals

  script:
  input = mdBam.collect{"-I $it"}.join(' ')

  """
  echo -e "idPatient:\t"${idPatient}"\nidSample:\t"${idSample}"\nmdBam:\t"${mdBam}"\n" > logInfo
  java -Xmx${task.memory.toGiga()}g -jar ${params.gatkHome}/GenomeAnalysisTK.jar \
  -T RealignerTargetCreator \
  $input \
  -R $gf \
  -known $ki \
  -known $mi \
  -nt ${task.cpus} \
  -o ${idPatient}.intervals
  """
}

intervals = logChannelContent("Intervals passed to Realign: ",intervals)

/*
 * realign, use nWayOut to split into tumor/normal again
 */

process Realign {

  memory { 16.GB * task.attempt }
  time { 20.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
  maxRetries 3
  maxErrors '-1'

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
  val(idPatient) into idPatient
  val(idSample) into tempSamples
  file("*.md.real.bam") into tempBams
  file("*.md.real.bai") into tempBais

  script:
  input = mdBam.collect{"-I $it"}.join(' ')

  """
  echo -e "idPatient:\t"${idPatient}"\nidSample:\t"${idSample}"\nmdBam:\t"${mdBam}"\nmdBai:\t"${mdBai}"\n" > logInfo
  java -Xmx${task.memory.toGiga()}g -jar ${params.gatkHome}/GenomeAnalysisTK.jar \
  -T IndelRealigner \
  $input \
  -R $gf \
  -targetIntervals $intervals \
  -known $ki \
  -known $mi \
  -nWayOut '.real.bam'
  """
}

// [maxime] If I make a set out of this process I got a list of lists, which cannot be iterate via a single process
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
realignedBam = idPatient.spread(tempSamples)

realignedBam = logChannelContent("realignedBam to BaseRecalibrator: ", realignedBam)

process CreateRecalibrationTable {

  module 'java/sun_jdk1.8.0_92'

  cpus 8
  memory { 8.GB * task.attempt }       // 6G is certainly low even for downsampled (30G) data
  time { 16.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
  maxRetries 3
  maxErrors '-1'

  input:
  set idPatient, idSample, file(realignedBamFile), file(realignedBaiFile) from realignedBam
  file refs["genomeFile"]
  file refs["dbsnp"]
  file refs["kgIndels"]
  file refs["millsIndels"]

  output:
  set idPatient, idSample, file(realignedBamFile), file(realignedBaiFile), file("${idSample}.recal.table") into recalibrationTable

  """
  java -Xmx${task.memory.toGiga()}g -Djava.io.tmpdir="/tmp" \
  -jar ${params.gatkHome}/GenomeAnalysisTK.jar \
  -T BaseRecalibrator \
  -R ${refs["genomeFile"]} \
  -I $realignedBamFile \
  -knownSites ${refs["dbsnp"]} \
  -knownSites ${refs["kgIndels"]} \
  -knownSites ${refs["millsIndels"]} \
  -nct ${task.cpus} \
  -l INFO \
  -o ${idSample}.recal.table
  """
}

recalibrationTable = logChannelContent("Base recalibrated table for recalibration: ", recalibrationTable)

process RecalibrateBam {

  memory { 16.GB * task.attempt }
  time { 16.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
  maxRetries 3
  maxErrors '-1'
  cpus 8

  input:
  set idPatient, idSample, file(realignedBamFile), file(realignedBaiFile), recalibrationReport from recalibrationTable
  file refs["genomeFile"]

  output:
  set idPatient, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bai") into recalibratedBams

  // TODO: ditto as at the previous BaseRecalibrator step, consider using -nct 4 
  """
  java -Xmx${task.memory.toGiga()}g -jar ${params.gatkHome}/GenomeAnalysisTK.jar \
  -T PrintReads \
  -R ${refs["genomeFile"]} \
  -nct ${task.cpus} \
  -I $realignedBamFile \
  --BQSR $recalibrationReport \
  -o ${idSample}.recal.bam
  """
}

recalibratedBams = logChannelContent("Recalibrated Bam for variant Calling: ", recalibratedBams)

// [maxime] Here we have a recalibrated bam set, but we need to separate the bam files based on patient status.
// The sample tsv config file which is formatted like: "subject status sample lane fastq1 fastq2"
// cf fastqFiles channel, I decided just to add __status to the sample name to have less changes to do.
// And so I'm sorting the channel if the sample match __0, then it's a normal sample, otherwise tumor.
// Then spread normal over tumor to get each possibilities
// ie. normal vs tumor1, normal vs tumor2, normal vs tumor3
// then copy this channel into channels for each variant calling
// I guess it will still work even if we have multiple normal samples

bamsTumor = Channel.create()
bamsNormal = Channel.create()

// separate recalibrate files by filename suffix: __0 means normal, __1 means tumor recalibrated BAM
recalibratedBams
  .choice(bamsTumor, bamsNormal) { it[1] =~ /__0$/ ? 1 : 0 }

bamsTumor = logChannelContent("Tumor Bam for variant Calling: ", bamsTumor)
bamsNormal = logChannelContent("Normal Bam for variant Calling: ", bamsNormal)

bamsAll = Channel.create()
bamsAll = bamsNormal.spread(bamsTumor)

// [maxime] Since idPatientNormal and idPatientTumor are the same, I'm removing it from BamsAll Channel
// I don't think a groupTuple can be used to do that, but it could be a good idea to look if there is a nicer way to do that

bamsAll = bamsAll.map {
  idPatientNormal, idSampleNormal, bamNormal, baiNormal, idPatientTumor, idSampleTumor, bamTumor, baiTumor ->
  [idPatientNormal, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor]
}

bamsAll = logChannelContent("Mapped Recalibrated Bam for variant Calling: ", bamsAll)

// [Szilva] We know that MuTect2 (and other somatic callers) are notoriously slow. To speed them up we are chopping the reference into 
// smaller pieces at centromeres (see repeates/centromeres.list), do variant calling by this intervals, and re-merge the VCFs.
// Since we are on a cluster, this can parallelize the variant call process, and push down the MuTect2 waiting time significanlty (1/10)

// first create channels for each variant caller
bamsForMuTect2 = Channel.create()
bamsForVarDict = Channel.create()

Channel
  .from bamsAll
  .separate( bamsForMuTect2, bamsForVarDict) {a -> [a, a]}

// define intervals file by --intervals
// TODO: add as a parameter file
intervalsFile = file(params.intervals)
intervals = Channel
    .from(intervalsFile.readLines())

// in fact we need two channels: one for the actual genomic region, and an other for names
// without ":", as nextflow is not happy with them (will report as a failed process).
// For region 1:1-2000 the output file name will be something like 1_1-2000_Sample_name.mutect2.vcf
// from the "1:1-2000" string make ["1:1-2000","1_1-2000"]
gI = intervals
  .map {a -> [a,a.replaceFirst(/\:/,"_")]}

MuTect2Intervals = Channel.create()
VarDictIntervals = Channel.create()
Channel
  .from gI
  .separate (MuTect2Intervals, VarDictIntervals) {a -> [a,a]}

// now add genomic intervals to the sample information
// join [idPatientNormal, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor] and ["1:1-2000","1_1-2000"] 
// and make a line for each interval

bamsFMT2 = bamsForMuTect2.spread(MuTect2Intervals)

process RunMutect2 {

  module 'bioinfo-tools'
  module 'java/sun_jdk1.8.0_92'

  threads 16
  memory { 16.GB * task.attempt }
  time { 16.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
  maxRetries 3
  maxErrors '-1'

  input:
  set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), genInt, gen_int from bamsFMT2

  output:
  set idPatient, val("${gen_int}_${idSampleNormal}_${idSampleTumor}"), file("${gen_int}_${idSampleNormal}_${idSampleTumor}.mutect2.vcf") into mutectVariantCallingOutput
  
  // we are using MuTect2 shipped in GATK v3.6
  """
  java -Xmx${task.memory.toGiga()}g -jar ${params.mutect2Home}/GenomeAnalysisTK.jar \
  -T MuTect2 \
  -nct ${task.threads} \
  -R ${refs["genomeFile"]} \
  --cosmic ${refs["cosmic"]} \
  --dbsnp ${refs["dbsnp"]} \
  -I:normal $bamNormal \
  -I:tumor $bamTumor \
  -L \"${genInt}\" \
  -o ${gen_int}_${idSampleNormal}_${idSampleTumor}.mutect2.vcf
  """
}

mutectVariantCallingOutput = logChannelContent("Mutect2 output: ", mutectVariantCallingOutput)
// TODO: merge call output

// we are doing the same trick for VarDictJava: running for the whole reference is a PITA, so we are chopping at repeats
// (or centromeres) where no useful variant calls are expected

bamsFVD = bamsForVarDict.spread(VarDictIntervals)

process VarDict {

  // ~/dev/VarDictJava/build/install/VarDict/bin/VarDict -G /sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/human_g1k_v37_decoy.fasta -f 0.1 -N "tiny" -b "tiny.tumor__1.recal.bam|tiny.normal__0.recal.bam" -z 1 -F 0x500 -c 1 -S 2 -E 3 -g 4 -R "1:131941-141339"
  // we need further filters, but some of the outputs are empty files, confusing the VCF generator script

  module 'bioinfo-tools'
  module 'VarDictJava/1.4.5'

  cpus 1
  memory { 16.GB * task.attempt }
  time { 16.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
  maxRetries 3
  maxErrors '-1'

  input:
  set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), genInt, gen_int from bamsFVD

  output:
  set idPatient, idSampleNormal, idSampleTumor, val("${gen_int}_${idSampleNormal}_${idSampleTumor}"), file("${gen_int}_${idSampleNormal}_${idSampleTumor}.VarDict.out") into varDictVariantCallingOutput

  """
  VarDict -G ${refs["genomeFile"]} \
  -f 0.01 -N $bamTumor \
  -b "$bamTumor|$bamNormal" \
  -z 1 -F 0x500 \
  -c 1 -S 2 -E 3 -g 4 \
  -R ${genInt} > ${gen_int}_${idSampleNormal}_${idSampleTumor}.VarDict.out
  """
}

// now we want to collate all the pieces of the VarDict outputs and concatenate the output files
// so we can feed them into the somatic filter and the VCF converter

varDictVariantCallingOutput = logChannelContent("VarDict VCF channel: ",varDictVariantCallingOutput)
(varDictVariantCallingOutput ,idPatient, idNormal, idTumor) = getPatientAndSample(varDictVariantCallingOutput)

vdFilePrefix = idPatient + "_" + idNormal + "_" + idTumor
vdFilesOnly = varDictVariantCallingOutput.map { x -> x.last()}

process VarDictCollatedVCF {
  publishDir "/home/szilva/dev/forkCAW/"

  module 'bioinfo-tools'
  module 'VarDictJava/1.4.5'
  module 'samtools/1.3'

  cpus 1
  memory { 16.GB * task.attempt }
  time { 16.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
  maxRetries 3
  maxErrors '-1'

  input:
  file vdPart from vdFilesOnly.toList()

  output:
  file(vdFilePrefix + ".VarDict.vcf")

  script:
  """
  for vdoutput in ${vdPart}
  do
    echo
    cat \$vdoutput | ${params.vardictHome}/testsomatic.R >> testsomatic.out
  done
  ${params.vardictHome}/var2vcf_somatic.pl -f 0.01 -N "${vdFilePrefix}" testsomatic.out > ${vdFilePrefix}.VarDict.vcf
  """
}

/*
========================================================================================
=                                   F U N C T I O N S                                  =
========================================================================================
*/

/* 
 * Helper function, given a file Path 
 * returns the file name region matching a specified glob pattern
 * starting from the beginning of the name up to last matching group.
 * 
 * For example: 
 *   readPrefix('/some/data/file_alpha_1.fa', 'file*_1.fa' )
 * 
 * Returns: 
 *   'file_alpha'
 */

def readPrefix (Path actual, template) {
  final fileName = actual.getFileName().toString()
  def filePattern = template.toString()
  int p = filePattern.lastIndexOf('/')
  if( p != -1 ) filePattern = filePattern.substring( p + 1 )
  if( !filePattern.contains('*') && !filePattern.contains('?') )
  filePattern = '*' + filePattern
  def regex = filePattern
    .replace('.', '\\.')
    .replace('*', '(.*)')
    .replace('?', '(.?)')
    .replace('{', '(?:')
    .replace('}', ')')
    .replace(',', '|')

  def matcher = (fileName =~ /$regex/)
  if ( matcher.matches() ) {
    def end = matcher.end( matcher.groupCount() )
    def prefix = fileName.substring(0,end)
    while ( prefix.endsWith('-') || prefix.endsWith('_') || prefix.endsWith('.') )
      prefix = prefix[0..-2]
    return prefix
  }
  return fileName
}

/*
 * Paolo Di Tommaso told that it should be solved by the Channel view() method, but frankly that was not working
 * as expected. This simple function makes two channels from one, and logs the content for one (thus consuming that channel)
 * and returns with the intact copy for further processing.
 */

def logChannelContent (aMessage, aChannel) {
  resChannel = Channel.create()
  logChannel = Channel.create()
  Channel
    .from aChannel
    .separate(resChannel,logChannel) {a -> [a, a]}
  logChannel.subscribe {log.info aMessage + " -- $it"}
  return resChannel
}

def getPatientAndSample(aCh) {
  consCh = Channel.create()
  originalCh = Channel.create()

  // get the patient ID
  // duplicate channel to get sample name
  Channel.from aCh.separate(consCh,originalCh) {x -> [x,x]}

  // use the "consumed" channel to get it
  // we are assuming the first column is the same for the patient, as hoping
  // people do not want to compare samples from different patients
  idPatient = consCh.map {x -> [x.get(0)]}.unique().getVal()[0]
  // we have to close to make sure remainding items are not
  consCh.close()

  // similar procedure for the normal sample name
  Channel.from originalCh.separate(consCh,originalCh) {x -> [x,x]}
  idNormal = consCh.map {x -> [x.get(1)]}.unique().getVal()[0]
  consCh.close()

  // ditto for the tumor
  Channel.from originalCh.separate(consCh,originalCh) {x -> [x,x]}
  idTumor = consCh.map {x -> [x.get(2)]}.unique().getVal()[0]
  consCh.close()

  return [originalCh, idPatient, idNormal, idTumor]
}

