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

String version = "0.0.3"
String dateUpdate = "2016-08-25"

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
 * Added a steps possibilities to do some process and skip others
 * borrowed the idea from https://github.com/guigolab/grape-nf
 */

switch (params) {
  case {params.help} :
    text = Channel.from(
      "CANCER ANALYSIS WORKFLOW ~ version $version",
      "    Usage:",
      "       nextflow run MultiFQtoVC.nf -c <file.config> --sample <sample.tsv>",
      "   [--steps STEP[,STEP]]",
        "       optional option for now, help you configure",
      "       which process will be processed by the workflow.",
      "       Possible values are:",
      "         preprocessing (default, will start workflow with FASTQ files)",
      "         nopreprocessing (will start workflow with recalibrated BAM files):",
      "         MuTect2 (use MuTect2 for VC):",
      "         VarDict (use VarDict for VC):",
      "         Strelka (use Strelka for VC):",
      "         Manta (use Manta for SV):",
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
  "genomeFile":   params.genome,       // genome reference
  "genomeIndex":  params.genomeIndex,  // genome reference index
  "genomeDict":   params.genomeDict,   // genome reference dictionary
  "kgIndels":     params.kgIndels,     // 1000 Genomes SNPs
  "kgIndex":      params.kgIndex,      // 1000 Genomes SNPs index
  "dbsnp":        params.dbsnp,        // dbSNP
  "dbsnpIndex":   params.dbsnpIndex,   // dbSNP index
  "millsIndels":  params.millsIndels,  // Mill's Golden set of SNPs
  "millsIndex":   params.millsIndex,   // Mill's Golden set index
  "sample":       params.sample,       // the sample sheet (multilane data refrence table, see below)
  "cosmic":       params.cosmic,       // cosmic vcf file
  "intervals":    params.intervals,    // intervals file for spread-and-gather processes (usually chromosome chunks at centromeres)
  "MantaRef":     params.mantaRef,     // copy of the genome reference file 
  "MantaIndex":   params.mantaIndex   // reference index indexed with samtools/0.1.19
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
 * Getting list of steps from comma-separated strings
 */

if (params.steps) {
  workflowSteps = params.steps.split(',').collect { it.trim() }
} else {
  workflowSteps = []
}

if ('preprocessing' in workflowSteps && 'nopreprocessing' in workflowSteps) {
  text = Channel.from(
    "CANCER ANALYSIS WORKFLOW ~ version $version",
    "Cannot use preprocessing and nopreprocessing at the same time")
  text.subscribe { println "$it" }
  exit 1
}

if (!('preprocessing' in workflowSteps && 'nopreprocessing' in workflowSteps)) {
  workflowSteps.add('preprocessing')
}

/*
 * Verifying the existence of the sample file.
 */

if (!params.sample) {
  text = Channel.from(
    "CANCER ANALYSIS WORKFLOW ~ version $version",
    "Missing the sample TSV config file: please specify it.",
    "    Usage",
    "       nextflow run MultiFQtoVC.nf -c <file.config> --sample <sample.tsv>")
  text.subscribe { println "$it" }
  exit 1
}

sampleTSVconfig = file(params.sample)

if ('preprocessing' in workflowSteps) {
  /*
   * Channeling the TSV file containing FASTQ for preprocessing
   * The format is: "subject status sample lane fastq1 fastq2"
   * ie: HCC1954 0 HCC1954.blood HCC1954.blood_1 data/HCC1954.normal_S22_L001_R1_001.fastq.gz data/HCC1954.normal_S22_L001_R2_001.fastq.gz
   * I just added __status to the idSample so that the whole pipeline is still working without having to change anything.
   * I know, it is lazy...
   */

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
  fastqFiles = logChannelContent("FASTQ files and IDs to process: ", fastqFiles)
} else {
  /*
   * Channeling the TSV file containing BAM for Recalibration
   * The format is: "subject status sample bam bai recal"
   * ie: HCC1954 0 HCC1954.blood data/HCC1954.normal.bam HCC1954.normal.bai HCC1954.normal.recal.table
   * Still with the __status added to the idSample
   */

  bamFiles = Channel
    .from(sampleTSVconfig.readLines())
    .map {line ->
      list        = line.split()
      idPatient   = list[0]
      idSample    = "${list[2]}__${list[1]}"
      bamFile     = file(list[3])
      baiFile     = file(list[4])
      recalTable  = file(list[5])
      [ idPatient, idSample, bamFile, baiFile, recalTable ]
  }
  bamFiles = logChannelContent("Bam files and IDs to process: ", bamFiles)
}

/*
========================================================================================
=                                   P R O C E S S E S                                  =
========================================================================================
*/

if ('preprocessing' in workflowSteps) {
  process Mapping {
    publishDir "Preprocessing/Mapping"

    module 'bioinfo-tools'
    module 'bwa/0.7.8'
    module 'samtools/1.3'

    cpus 8
    memory { 16.GB * task.attempt }
    time { 20.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

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

  bams = logChannelContent("BAM files before sorting into group or single:", bams)

  /*
   * Borrowed code from https://github.com/guigolab/chip-nf
   * Now, we decide whether bam is standalone or should be merged by sample (id (column 1) from channel bams)
   * http://www.nextflow.io/docs/latest/operator.html?highlight=grouptuple#grouptuple
   */

  /*
   * Merge or rename bam
   */

  singleBam = Channel.create()
  groupedBam = Channel.create()

  bams.groupTuple(by:[1])
    .choice(singleBam, groupedBam) { it[3].size() > 1 ? 1 : 0 }

  singleBam = logChannelContent("Single BAMs before merge:", singleBam)
  groupedBam = logChannelContent("Grouped BAMs before merge:", groupedBam)

  process MergeBam {
    publishDir "Preprocessing/MergeBam"

    module 'bioinfo-tools'
    module 'picard/1.118'

    memory { 16.GB * task.attempt }
    time { 16.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    input:
    set idPatient, idSample, idRun, file(bam) from groupedBam

    output:
    set idPatient, idSample, idRun, file("${idSample}.bam") into mergedBam

    script:
		// make a string list, replace space and comma as ".bam" and feed the string to shell
		// so from a run ID string like
		// "[tiny.normal_4, tiny.normal_2, tiny.normal_7, tiny.normal_1, tiny.normal_8]"
		// we will have
		// "tiny.normal_4.bam INPUT=tiny.normal_2.bam INPUT=tiny.normal_7.bam INPUT=tiny.normal_1.bam INPUT=tiny.normal_8.bam"
		bamInputListStr = idRun.toListString()
		bamInput = bamInputListStr.replace(', ', '.bam INPUT=').replace(']','.bam').replace('[','')
    idRun = idRun.sort().join(':')
    """
    echo -e "idPatient:\t"${idPatient}"\nidSample:\t"${idSample}"\nidRun:\t"${idRun}"\nbam:\t"${bam}"\n" > logInfo
		BAM_INPUT=${bamInput}
    java -Xmx${task.memory.toGiga()}g -jar ${params.picardHome}/MergeSamFiles.jar \
    INPUT=${bamInput} \
    TMP_DIR=. \
    ASSUME_SORTED=true \
		USE_THREADING=true \
    VALIDATION_STRINGENCY=LENIENT \
    CREATE_INDEX=FALSE \
    OUTPUT=${idSample}.bam
    """
  }

  // Renaming is totally useless, but the file name is consistent with the rest of the pipeline

  process RenameSingleBam {
    publishDir "Preprocessing/RenameSingleBam"

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
    publishDir "Preprocessing/MarkDuplicates"

    module 'bioinfo-tools'
    module 'picard/1.118'

    memory { 16.GB * task.attempt }
    time { 16.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

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
  duplicatesRealign = Channel.create()
  duplicatesRealign = duplicatesForRealignement.groupTuple()
  duplicatesRealign = logChannelContent("BAMs for IndelRealigner grouped by overall subject/patient ID: ", duplicatesRealign)

  /*
   * Creating target intervals for indel realigner.
   * Though VCF indexes are not needed explicitly, we are adding them so they will be linked, and not re-created on the fly.
   */

  process CreateIntervals {
    publishDir "Preprocessing/CreateIntervals"

    module 'java/sun_jdk1.8.0_92'

    cpus 8
    memory { 16.GB * task.attempt }
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
    publishDir "Preprocessing/Realign"

    module 'java/sun_jdk1.8.0_92'

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
  realignedBam = idPatient.spread(tempSamples)

  realignedBam = logChannelContent("realignedBam to BaseRecalibrator: ", realignedBam)

  process CreateRecalibrationTable {
    publishDir "Preprocessing/CreateRecalibrationTable"

    module 'java/sun_jdk1.8.0_92'

    cpus 8
    memory { 16.GB * task.attempt }
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
} else {
  recalibrationTable = bamFiles
}

process RecalibrateBam {
  publishDir "Preprocessing/RecalibrateBam"

  module 'java/sun_jdk1.8.0_92'

  cpus 8
  memory { 16.GB * task.attempt }
  time { 16.h * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
  maxRetries 3
  maxErrors '-1'

  input:
  set idPatient, idSample, file(realignedBamFile), file(realignedBaiFile), recalibrationReport from recalibrationTable
  file refs["genomeFile"]

  output:
  set idPatient, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bai") into recalibratedBams

  // TODO: ditto as at the previous BaseRecalibrator step, consider using -nct 4 
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar ${params.gatkHome}/GenomeAnalysisTK.jar \
  -T PrintReads \
  -R ${refs["genomeFile"]} \
  -nct ${task.cpus} \
  -I $realignedBamFile \
  --BQSR $recalibrationReport \
  -o ${idSample}.recal.bam
  """
}

recalibratedBams = logChannelContent("Recalibrated Bam for variant Calling: ", recalibratedBams)

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

// separate recalibrate files by filename suffix: __0 means normal, __1 means tumor recalibrated BAM
recalibratedBams
  .choice(bamsTumor, bamsNormal) { it[1] =~ /__0$/ ? 1 : 0 }

bamsTumor = logChannelContent("Tumor Bam for variant Calling: ", bamsTumor)
bamsNormal = logChannelContent("Normal Bam for variant Calling: ", bamsNormal)

bamsAll = Channel.create()
bamsAll = bamsNormal.spread(bamsTumor)

// Since idPatientNormal and idPatientTumor are the same, I'm removing it from BamsAll Channel
// I don't think a groupTuple can be used to do that, but it could be a good idea to look if there is a nicer way to do that

bamsAll = bamsAll.map {
  idPatientNormal, idSampleNormal, bamNormal, baiNormal, idPatientTumor, idSampleTumor, bamTumor, baiTumor ->
  [idPatientNormal, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor]
}

bamsAll = logChannelContent("Mapped Recalibrated Bam for variant Calling: ", bamsAll)

// We know that MuTect2 (and other somatic callers) are notoriously slow. To speed them up we are chopping the reference into 
// smaller pieces at centromeres (see repeates/centromeres.list), do variant calling by this intervals, and re-merge the VCFs.
// Since we are on a cluster, this can parallelize the variant call process, and push down the variant call wall clock time significanlty.

// first create channels for each variant caller
bamsForMuTect2 = Channel.create()
bamsForVarDict = Channel.create()
bamsForManta = Channel.create()
bamsForStrelka = Channel.create()

Channel
  .from bamsAll
  .separate(bamsForMuTect2, bamsForVarDict, bamsForManta, bamsForStrelka) {a -> [a, a, a, a]}

// define intervals file by --intervals
intervalsFile = file(params.intervals)

intervals = Channel
  .from(intervalsFile.readLines())

// in fact we need two channels: one for the actual genomic region, and an other for names
// without ":", as nextflow is not happy with them (will report as a failed process).
// For region 1:1-2000 the output file name will be something like 1_1-2000_Sample_name.mutect2.vcf
// from the "1:1-2000" string make ["1:1-2000","1_1-2000"]
gI = intervals
  .map {a -> [a,a.replaceFirst(/\:/,"_")]}

muTect2Intervals = Channel.create()
varDictIntervals = Channel.create()
strelkaIntervals = Channel.create()

Channel
  .from gI
  .separate (muTect2Intervals, varDictIntervals, strelkaIntervals) {a -> [a, a, a]}

// now add genomic intervals to the sample information
// join [idPatientNormal, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor] and ["1:1-2000","1_1-2000"] 
// and make a line for each interval

if ('MuTect2' in workflowSteps) {

  bamsFMT2 = bamsForMuTect2.spread(muTect2Intervals)
  bamsFMT2 = logChannelContent("Bams for Mutect2: ", bamsFMT2)

  process RunMutect2 {
    publishDir "VariantCalling/MuTect2/intervals"

    module 'bioinfo-tools'
    module 'java/sun_jdk1.8.0_92'

    cpus 8 
    memory { 16.GB * task.attempt }
    time { 16.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    input:
    set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), genInt, gen_int from bamsFMT2

    output:
    set idPatient, idSampleNormal, idSampleTumor, val("${gen_int}_${idSampleNormal}_${idSampleTumor}"), file("${gen_int}_${idSampleNormal}_${idSampleTumor}.mutect2.vcf") into mutect2VariantCallingOutput
  
    // we are using MuTect2 shipped in GATK v3.6
    """
    java -Xmx${task.memory.toGiga()}g -jar ${params.mutect2Home}/GenomeAnalysisTK.jar \
    -T MuTect2 \
    -nct ${task.cpus} \
    -R ${refs["genomeFile"]} \
    --cosmic ${refs["cosmic"]} \
    --dbsnp ${refs["dbsnp"]} \
    -I:normal $bamNormal \
    -I:tumor $bamTumor \
    -L \"${genInt}\" \
    -o ${gen_int}_${idSampleNormal}_${idSampleTumor}.mutect2.vcf
    """
  }

  // we are expecting one patient, one normal, and usually one, but occasionally more than one tumor
  // samples (i.e. relapses). The actual calls are always related to the normal, but spread across
  // different intervals. So, we have to collate (merge) intervals for each tumor case if there are
  // more than one. Therefore, what we want to do is to filter the multiple tumor cases into separate 
  // channels and collate them according to their stage.
  mutect2VariantCallingOutput = logChannelContent("Mutect2 output: ", mutect2VariantCallingOutput)
  filesToCollate = mutect2VariantCallingOutput
  .groupTuple(by: 2)
  .map { 
    x ->  [
      x[0].get(0),  // the patient ID
      x[1].get(0),  // ID of the normal sample 
      x[2],         // ID of the tumor sample (primary, relapse, whatever)
      x[4]          // list of VCF files
      ]
    }

  // we have to separate IDs and files
  collatedIDs = Channel.create()
  collatedFiles = Channel.create()
  tumorEntries = Channel.create()
  Channel
    .from filesToCollate
    .separate(collatedIDs, collatedFiles, tumorEntries) {x -> [ x, [x[0],x[1],x[2]], x[2] ]}

  (idPatient, idNormal) = getPatientAndNormalIDs(collatedIDs)
  println "Patient's ID: " + idPatient
  println "Normal ID: " + idNormal
  pd = "VariantCalling/MuTect2"
  process concatFiles {
    publishDir = pd

    module 'bioinfo-tools'
    module 'java/sun_jdk1.8.0_92'

    cpus 8 
    memory { 16.GB * task.attempt }
    time { 16.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    input:
    set idT from tumorEntries

    output:
    file "MuTect2*.vcf"

    script:
    """
    VARIANTS=`ls ${workflow.launchDir}/${pd}/intervals/*${idT}*.mutect2.vcf| awk '{printf(" -V %s\\n",\$1) }'`
    java -Xmx${task.memory.toGiga()}g -cp ${params.mutect2Home}/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R ${refs["genomeFile"]}  \$VARIANTS -out MuTect2_${idPatient}_${idNormal}_${idT}.vcf
    """
  }
}

if ('VarDict' in workflowSteps) {
  // we are doing the same trick for VarDictJava: running for the whole reference is a PITA, so we are chopping at repeats
  // (or centromeres) where no useful variant calls are expected
  bamsFVD = bamsForVarDict.spread(varDictIntervals)
  bamsFVD = logChannelContent("Bams for VarDict: ", bamsFVD)

  process VarDict {
    publishDir "VariantCalling/VarDictJava"
  
    // ~/dev/VarDictJava/build/install/VarDict/bin/VarDict -G /sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/human_g1k_v37_decoy.fasta -f 0.1 -N "tiny" -b "tiny.tumor__1.recal.bam|tiny.normal__0.recal.bam" -z 1 -F 0x500 -c 1 -S 2 -E 3 -g 4 -R "1:131941-141339"
    // we need further filters, but some of the outputs are empty files, confusing the VCF generator script
  
    module 'bioinfo-tools'
    module 'java/sun_jdk1.8.0_92'
    module 'R/3.2.3'
    module 'gcc/4.9.2'
    module 'java/sun_jdk1.8.0_40'
    module 'perl/5.18.4'

    cpus 1
    memory { 6.GB * task.attempt }
    time { 16.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'
  
    input:
    set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), genInt, gen_int from bamsFVD
  
    output:
    set idPatient, idSampleNormal, idSampleTumor, val("${gen_int}_${idSampleNormal}_${idSampleTumor}"), file("${gen_int}_${idSampleNormal}_${idSampleTumor}.VarDict.out") into varDictVariantCallingOutput

    """
    ${params.vardictHome}/vardict.pl -G ${refs["genomeFile"]} \
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
  
  resultsDir = file("${idPatient}.results")
  resultsDir.mkdir() 
  
  process VarDictCollatedVCF {
    publishDir "VariantCalling/VarDictJava"
  
    module 'bioinfo-tools'
    module 'samtools/1.3'
    module 'java/sun_jdk1.8.0_92'
    module 'R/3.2.3'
    module 'gcc/4.9.2'
    module 'java/sun_jdk1.8.0_40'
    module 'perl/5.18.4'
  
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
} else {
  bamsForVarDict.close()
  varDictIntervals.close()
}

if ('Strelka' in workflowSteps) {
  bamsForStrelka = logChannelContent("Bams for Strelka: ", bamsForStrelka)
  strelkaIntervals = logChannelContent("Intervals for Strelka: ", strelkaIntervals)
  bamsFSTR = bamsForStrelka.spread(strelkaIntervals)
  bamsFSTR = logChannelContent("Bams with Intervals for Strelka: ", bamsFSTR)

  process RunStrelka {
    publishDir "VariantCalling/Strelka"

    module 'bioinfo-tools'

    cpus 1
    memory { 16.GB * task.attempt }
    time { 16.h * task.attempt }
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
    maxRetries 3
    maxErrors '-1'

    input:
    set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), genInt, gen_int from bamsFSTR
    file refs["genomeFile"]
    file refs["genomeIndex"]

    output:
    sed idPatient, val("${gen_int}_${idSampleNormal}_${idSampleTumor}"), file("*.vcf") into strelkaVariantCallingOutput

    """
    ${params.strelkaHome}/bin/configureStrelkaWorkflow.pl \
    --tumor ${bamTumor} \
    --normal ${bamNormal} \
    --ref ${refs["genomeFile"]} \
    --config ${params.strelkaCFG} \
    --output-dir strelka

    cd strelka

    make -j 16
    """
  }
  strelkaVariantCallingOutput = logChannelContent("Strelka output: ", strelkaVariantCallingOutput)
} else {
  bamsForStrelka = logChannelContent("Bams for Strelka: ", bamsForStrelka)
  bamsForStrelka.close()
  strelkaIntervals.close()
}

if ('Manta' in workflowSteps) {
  process Manta {
    module 'bioinfo-tools'
    module 'manta/0.27.1'
    module 'samtools/0.1.19'

    cpus 8

    input:
        file refs["MantaRef"]
        file refs["MantaIndex"]
        set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from bamsForManta
    
    output:
       set idPatient, val("${idSampleNormal}_${idSampleTumor}"),file("${idSampleNormal}_${idSampleTumor}.somaticSV.vcf"),file("${idSampleNormal}_${idSampleTumor}.candidateSV.vcf"),file("${idSampleNormal}_${idSampleTumor}.diploidSV.vcf"),file("${idSampleNormal}_${idSampleTumor}.candidateSmallIndels.vcf")  into mantaVariantCallingOutput


    //NOTE: Manta is very picky about naming and reference indexes, the input bam should not contain too many _ and the reference index must be generated using a supported samtools version.
    //Moreover, the bam index must be named .bam.bai, otherwise it will not be recognized

    """
    mv ${bamNormal} Normal.bam
    mv ${bamTumor} Tumor.bam

    mv ${baiNormal} Normal.bam.bai
    mv ${baiTumor} Tumor.bam.bai

    configManta.py --normalBam Normal.bam --tumorBam Tumor.bam --reference ${refs["MantaRef"]} --runDir MantaDir
    python MantaDir/runWorkflow.py -m local -j 8
    gunzip -c MantaDir/results/variants/somaticSV.vcf.gz > ${idSampleNormal}_${idSampleTumor}.somaticSV.vcf
    gunzip -c MantaDir/results/variants/candidateSV.vcf.gz > ${idSampleNormal}_${idSampleTumor}.candidateSV.vcf
    gunzip -c MantaDir/results/variants/diploidSV.vcf.gz > ${idSampleNormal}_${idSampleTumor}.diploidSV.vcf
    gunzip -c MantaDir/results/variants/candidateSmallIndels.vcf.gz > ${idSampleNormal}_${idSampleTumor}.candidateSmallIndels.vcf
    """
  }
  mantaVariantCallingOutput = logChannelContent("Manta output: ", mantaVariantCallingOutput)
} else {
  bamsForManta.close()
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

  aCh = logChannelContent("Channel content: ", aCh)
  patientsCh = Channel.create()
  normalCh = Channel.create()
  tumorCh = Channel.create()
  originalCh = Channel.create()

  // get the patient ID
  // duplicate channel to get sample name
  Channel.from aCh.separate(patientsCh, normalCh, tumorCh, originalCh) {x -> [x,x,x,x]} 

  // consume the patiensCh channel to get the patient's ID
  // we are assuming the first column is the same for the patient, as hoping 
  // people do not want to compare samples from differnet patients
  idPatient = patientsCh.map { x -> [x.get(0)]}.unique().getVal()[0]
  // we have to close to make sure remainding items are not waiting
  patientsCh.close()

  idNormal = normalCh.map { x -> [x.get(1)]}.unique().getVal()[0]  
  normalCh.close()

  idTumor = tumorCh.map { x -> [x.get(2)]}.unique().getVal()[0]  
  tumorCh.close()

  return [ originalCh, idPatient, idNormal, idTumor]
}

// TODO: merge the VarDict and the MuTect2 part
def getPatientAndNormalIDs(aCh) {
  patientsCh = Channel.create()
  normalCh = Channel.create()

  // get the patient ID
  // multiple the channel to get sample name
  Channel.from aCh.separate(patientsCh, normalCh) {x -> [x,x]} 

  // consume the patiensCh channel to get the patient's ID
  // we are assuming the first column is the same for the patient, as hoping 
  // people do not want to compare samples from differnet patients
  idPatient = patientsCh.map { x -> [x.get(0)]}.unique().getVal()[0]
  // we have to close to make sure remainding items are not waiting
  patientsCh.close()
  // something similar for the normal ID
  idNormal = normalCh.map { x -> [x.get(1)]}.unique().getVal()[0]  
  normalCh.close()

  return [idPatient, idNormal]
}
