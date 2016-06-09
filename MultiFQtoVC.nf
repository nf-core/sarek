#!/usr/bin/env nextflow

/*
 * Cancer Analysis Workflow for tumor/normal samples. Usage :
 *
 * $ nextflow run MultiFQtoVC.nf -c <file.config> --sample <sample.tsv>
 *
 * Parameters (like location of the genome reference, dbSNP location...) have to be given at the command-line.
 * We are not giving parameters here, but checking the existence of files.
 *
 * You also must give the sample.tsv file that is defining the tumor/normal pairs. See explanation below.
 */

//
// We want to list all the missing parameters at the beginning, hence checking only this value
// and exiting only at the very end if some of the parameters failed
//

// ############################### CONFIGURATION ###############################

String version    = "0.0.1"
String dateUpdate = "2016-06-08"

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
      "    Usage",
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

parametersDefined = true

/*
 * Use this closure to loop through all the parameters.
 * We can get an AssertionError exception from the file() method as well.
 */

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
  "genomeFile":     params.genome,      // genome reference
  "genomeIndex":    params.genomeIndex, // genome reference index
  "genomeDict":     params.genomeDict,  // genome reference dictionary
  "kgIndels":       params.kgIndels,    // 1000 Genomes SNPs
  "kgIndex":        params.kgIndex,     // 1000 Genomes SNPs index
  "dbsnp":          params.dbsnp,       // dbSNP
  "dbsnpIndex":     params.dbsnpIndex,  // dbSNP index
  "millsIndels":    params.millsIndels, // Mill's Golden set of SNPs
  "millsIndex":     params.millsIndex,  // Mill's Golden set index
  "sample":         params.sample       // the sample sheet (multilane data refrence table, see below)
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
 * Time to check the sample file. Its format is like: "subject sample lane fastq1 fastq2":
 * tcga.cl	tcga.cl.normal	tcga.cl.normal_1	data/tcga.cl.normal_L001_R1.fastq.gz	data/tcga.cl.normal_L001_R2.fastq.gz
 * tcga.cl	tcga.cl.tumor	tcga.cl.tumor_1	data/tcga.cl.tumor_L001_R1.fastq.gz	data/tcga.cl.tumor_L001_R2.fastq.gz
 * tcga.cl	tcga.cl.tumor	tcga.cl.tumor_2	data/tcga.cl.tumor_L002_R1.fastq.gz	data/tcga.cl.tumor_L002_R2.fastq.gz
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
 * Read config file, lets presume its "subject sample fastq1 fastq2"
 * for now and channel this out for mapping
 */

fastqFiles = Channel
  .from(sampleTSVconfig.readLines())
  .map { line ->
    list        = line.split()
    idPatient   = list[0]
    idSample    = list[1]
    idRun       = list[2]
    fastqFile1  = file(list[3])
    fastqFile2  = file(list[4])
    [ idPatient, idSample, idRun, fastqFile1, fastqFile2 ]
}

// ################################# PROCESSES #################################

fastqFiles = logChannelContent("FASTQ files and IDs to process: ",fastqFiles)

process Mapping {

  module 'bioinfo-tools'
  module 'bwa/0.7.8'
  module 'samtools/1.3'

  cpus 1

  input:
  file refs["genomeFile"]
  set idPatient, idSample, idRun, file(fq1), file(fq2) from fastqFiles

  output:
  set idPatient, idSample, idRun, file("${idRun}.bam") into bams

  // here I use params.genome for bwa ref so I dont have to link to all bwa index files

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
 * Borrowed code from chip.nf
 *
 * Now, we decide whether bam is standalone or should be merged by sample (id (column 1) from channel bams)
 * http://www.nextflow.io/docs/latest/operator.html?highlight=grouptuple#grouptuple
 *
 */

// Merge or rename bam
// Renaming is totally useless, but it is more consistent

singleBam  = Channel.create()
groupedBam = Channel.create()

bams.groupTuple(by:[1])
  .choice(singleBam, groupedBam) { it[3].size() > 1 ? 1 : 0 }

singleBam  = logChannelContent("Single BAMs before merge:", singleBam)
groupedBam = logChannelContent("Grouped BAMs before merge:", groupedBam)

process MergeBam {

  module 'bioinfo-tools'
  module 'samtools/1.3'

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
mergedBam        = logChannelContent("GROUPED: ", mergedBam)

/*
 * merge all bams (merged and singles) to a single channel
 */

bamList = Channel.create()
bamList = mergedBam.mix(singleRenamedBam)
// bamList = logChannelContent("Mixed channels: ", bamList)
bamList = bamList.map { idPatient, idSample, idRun, bam -> [idPatient[0], idSample, bam].flatten() }
// bamList = logChannelContent("Mapped and flattened: ", bamList)

bamList = logChannelContent("BAM list for MarkDuplicates: ",bamList)

/*
 *  mark duplicates all bams
 */

process MarkDuplicates {

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
  java -Xmx7g -jar ${params.picardHome}/MarkDuplicates.jar \
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

duplicatesForInterval = logChannelContent("BAMs for IndelRealigner before groupTuple: ",  duplicatesForInterval)

// group the marked duplicates Bams intervals by overall subject/patient id (idPatient)
duplicatesInterval = Channel.create()
duplicatesInterval = duplicatesForInterval.groupTuple()
duplicatesInterval = logChannelContent("BAMs for RealignerTargetCreator grouped by overall subject/patient ID: ",  duplicatesInterval)

duplicatesForRealignement = logChannelContent("BAMs for IndelRealigner before groupTuple: ",  duplicatesForRealignement)

// group the marked duplicates Bams for realign by overall subject/patient id (idPatient)
duplicatesRealign  = Channel.create()
duplicatesRealign  = duplicatesForRealignement.groupTuple()
duplicatesRealign  = logChannelContent("BAMs for IndelRealigner grouped by overall subject/patient ID: ",  duplicatesRealign)

/*
 * Creating target intervals for indel realigner.
 * Though VCF indexes are not needed explicitly, we are adding them so they will be linked, and not re-created on the fly.
 */

process CreateIntervals {

  cpus 6

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
  java -Xmx7g -jar ${params.gatkHome}/GenomeAnalysisTK.jar \
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
  set idPatient, idSample, file("*.md.real.bam"), file("*.md.real.bai") into realignedBam
  val(idPatient) into idPatient

  script:
  input = mdBam.collect{"-I $it"}.join(' ')

  """
  echo -e "idPatient:\t"${idPatient}"\nidSample:\t"${idSample}"\nmdBam:\t"${mdBam}"\nmdBai:\t"${mdBai}"\n" > logInfo
  java -Xmx7g -jar ${params.gatkHome}/GenomeAnalysisTK.jar \
  -T IndelRealigner \
  $input \
  -R $gf \
  -targetIntervals $intervals \
  -known $ki \
  -known $mi \
  -nWayOut '.real.bam'
  """
}

tempSamples  = Channel.create()
tempBams     = Channel.create()
tempBais     = Channel.create()

Channel
  .from realignedBam
  .separate( tempSamples, tempBams, tempBais ) { a -> [a, a, a] }

tempSamples  = tempSamples.flatMap { id, samples, bams, bais -> [samples] }.flatten().toSortedList().flatten()
tempBams     = tempBams.flatMap { id, samples, bams, bais -> [bams] }.flatten().toSortedList().flatten()
tempBais     = tempBais.flatMap { id, samples, bams, bais -> [bais] }.flatten().toSortedList().flatten()
tempSamples  = tempSamples.merge( tempBams, tempBais ) { s, b, i -> [s, b, i] }
realignedBam  = idPatient.spread(tempSamples)

realignedBam = logChannelContent("realignedBam to BaseRecalibrator: ", realignedBam)

process CreateRecalibrationTable {

  cpus 2

   input:
   set idPatient, idSample, realignedBamTable, realignedBaiTable from realignedBam
   file refs["genomeFile"]
   file refs["dbsnp"]
   file refs["kgIndels"]
   file refs["millsIndels"]

   output:
   set idPatient, idSample, realignedBamTable, file("${idSample}.recal.table") into recalibrationTable

   """
   java -Xmx7g -Djava.io.tmpdir=\$SNIC_TMP \
   -jar ${params.gatkHome}/GenomeAnalysisTK.jar \
   -T BaseRecalibrator \
   -R ${refs["genomeFile"]} \
   -I $realignedBamTable \
   -knownSites ${refs["dbsnp"]} \
   -knownSites ${refs["kgIndels"]} \
   -knownSites ${refs["millsIndels"]} \
   -nct ${task.cpus} \
   -l INFO \
   -o ${idSample}.recal.table
   """
}

recalibrationTable = logChannelContent("Base recalibrated table for recalibration: ",recalibrationTable)

process RecalibrateBam {

  input:
  set idPatient, idSample, recalibrationReport, recalibratedTable into recalibrationTable
  file refs["genomeFile"]
  file refs["dbsnp"]
  file refs["kgIndels"]
  file refs["millsIndels"]

  output:
  set idPatient, idSample, file("${idSample}.recal.bam") into recalibratedBam
  file '*.recal.bam' into recal_bam
  file '*.recal.bai' into recal_bai

  // java -jar GenomeAnalysisTK.jar \
  // -T PrintReads \
  // -R reference.fasta \
  // -I input.bam \
  // -BQSR recalibration_report.grp \
  // -o output.bam

  """
  java -Xmx7g -jar ${params.gatkHome}/GenomeAnalysisTK.jar \
  -T PrintReads \
  -R ${refs["genomeFile"]} \
  -I $recalibratedTable \
  --BQSR $recalibrationReport \
  -o ${idSample}.recal.bam
  """
}

recalibratedBam = logChannelContent("Recalibrated Bam for variant Calling: ",recalibratedBam)

/*
 * Variant Calling
 * The idea here is to make variant calls
 * normal/tumor
 * or if several tumors available
 * normal/tumor_1
 * normal/tumor_2
 * normal/tumor_3
 * ...
 *
 * IDEAS :
 * Duplicate output from Realign
 * and filter on name
 *
 */

// Here the difficulty will be to get the difference between normal/tumor
// cf https://groups.google.com/forum/m/#!topic/nextflow/zId8iRE6Z2w

// separate bams and inputs
// treat = Channel.create()
// control = Channel.create()
// bams.choice(treat, control) {
//   it[3] == 'input' ? 1 : 0
// }


// ################################# FUNCTIONS #################################

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
    .separate(resChannel,logChannel) { a -> [a, a] }
  logChannel.subscribe { log.info aMessage + " -- $it" }
  return resChannel
}