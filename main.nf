#!/usr/bin/env nextflow

/*
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
 - MarkDuplicates - Mark Duplicates with GATK4
 - CreateRecalibrationTable - Create Recalibration Table with BaseRecalibrator
 - RecalibrateBam - Recalibrate Bam with PrintReads
 - RunSamtoolsStats - Run Samtools stats on recalibrated BAM files
 - RunBamQCmapped - Run qualimap BamQC on mapped BAM files
 - RunBamQCrecalibrated - Run qualimap BamQC on recalibrated BAM files
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

// Check for awsbatch profile configuration
// make sure queue is defined
if (workflow.profile == 'awsbatch') {
    if(!params.awsqueue) exit 1, "Provide the job queue for aws batch!"
}


if (params.help) exit 0, helpMessage()
if (!SarekUtils.isAllowedParams(params)) exit 1, "params unknown, see --help for more information"
if (!checkUppmaxProject()) exit 1, "No UPPMAX project ID found! Use --project <UPPMAX Project ID>"

step = params.step.toLowerCase()
if (step == 'preprocessing') step = 'mapping'

directoryMap = SarekUtils.defineDirectoryMap(params.outDir)
referenceMap = defineReferenceMap()
stepList = defineStepList()

if (!SarekUtils.checkParameterExistence(step, stepList)) exit 1, 'Unknown step, see --help for more information'
if (step.contains(',')) exit 1, 'You can choose only one step, see --help for more information'
if (step == 'mapping' && !checkExactlyOne([params.test, params.sample, params.sampleDir]))
  exit 1, 'Please define which samples to work on by providing exactly one of the --test, --sample or --sampleDir options'
if (!SarekUtils.checkReferenceMap(referenceMap)) exit 1, 'Missing Reference file(s), see --help for more information'

if (params.test && params.genome in ['GRCh37', 'GRCh38']) {
  referenceMap.intervals = file("$workflow.projectDir/repeats/tiny_${params.genome}.list")
}

tsvPath = ''
if (params.sample) tsvPath = params.sample

// No need for tsv file for step annotate
if (!params.sample && !params.sampleDir) {
  tsvPaths = [
      'mapping': "${workflow.projectDir}/Sarek-data/testdata/tsv/tiny.tsv",
      'recalibrate': "${directoryMap.duplicateMarked}/duplicateMarked.tsv"
  ]
  if (params.test || step != 'mapping') tsvPath = tsvPaths[step]
}

// Set up the fastqFiles and bamFiles channels. One of them remains empty
fastqFiles = Channel.empty()
bamFiles = Channel.empty()
if (tsvPath) {
  tsvFile = file(tsvPath)
  switch (step) {
    case 'mapping': fastqFiles = extractFastq(tsvFile); break
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

if (step == 'mapping') (patientGenders, fastqFiles) = SarekUtils.extractGenders(fastqFiles)
else (patientGenders, bamFiles) = SarekUtils.extractGenders(bamFiles)

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

  publishDir "${directoryMap.fastQC}/${idRun}", mode: params.publishDirMode

  input:
    set idPatient, status, idSample, idRun, file(fastqFile1), file(fastqFile2) from fastqFilesforFastQC

  output:
    file "*_fastqc.{zip,html}" into fastQCreport

  when: step == 'mapping' && !params.noReports

  script:
  """
  fastqc -t 2 -q ${fastqFile1} ${fastqFile2}
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
    set idPatient, status, idSample, idRun, file("${idRun}.bam") into (mappedBam, mappedBamForQC)

  when: step == 'mapping' && !params.onlyQC

  script:
  CN = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ""
  readGroup = "@RG\\tID:${idRun}\\t${CN}PU:${idRun}\\tSM:${idSample}\\tLB:${idSample}\\tPL:illumina"
  // adjust mismatch penalty for tumor samples
  extra = status == 1 ? "-B 3" : ""
  """
  bwa mem -R \"${readGroup}\" ${extra} -t ${task.cpus} -M \
  ${genomeFile} ${fastqFile1} ${fastqFile2} | \
  samtools sort --threads ${task.cpus} -m 2G - > ${idRun}.bam
  """
}

if (params.verbose) mappedBam = mappedBam.view {
  "Mapped BAM (single or to be merged):\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\tRun   : ${it[3]}\n\
  File  : [${it[4].fileName}]"
}

process RunBamQCmapped {
  tag {idPatient + "-" + idSample}

  publishDir directoryMap.bamQC, mode: params.publishDirMode

  input:
    set idPatient, status, idSample, idRun, file(bam) from mappedBamForQC

  output:
    file(idSample) into bamQCmappedReport

  when: !params.noReports && !params.noBAMQC

  script:
  """
  qualimap --java-mem-size=${task.memory.toGiga()}G \
  bamqc \
  -bam ${bam} \
  --paint-chromosome-limits \
  --genome-gc-distr HUMAN \
  -nt ${task.cpus} \
  -skip-duplicated \
  --skip-dup-mode 0 \
  -outdir ${idSample} \
  -outformat HTML
  """
}

if (params.verbose) bamQCmappedReport = bamQCmappedReport.view {
  "BamQC report:\n\
  Dir   : [${it.fileName}]"
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

  publishDir params.outDir, mode: params.publishDirMode,
    saveAs: {
      if (it == "${idSample}.bam.metrics") "${directoryMap.markDuplicatesQC.minus(params.outDir+'/')}/${it}"
      else "${directoryMap.duplicateMarked.minus(params.outDir+'/')}/${it}"
    }

  input:
    set idPatient, status, idSample, file("${idSample}.bam") from mergedBam

  output:
    set idPatient, file("${idSample}_${status}.md.bam"), file("${idSample}_${status}.md.bai") into duplicateMarkedBams
    set idPatient, status, idSample, val("${idSample}_${status}.md.bam"), val("${idSample}_${status}.md.bai") into markDuplicatesTSV
    file ("${idSample}.bam.metrics") into markDuplicatesReport

  when: step == 'mapping' && !params.onlyQC

  script:
  """
  gatk --java-options ${params.markdup_java_options} \
  MarkDuplicates \
  --MAX_RECORDS_IN_RAM 50000 \
  --INPUT ${idSample}.bam \
  --METRICS_FILE ${idSample}.bam.metrics \
  --TMP_DIR . \
  --ASSUME_SORT_ORDER coordinate \
  --CREATE_INDEX true \
  --OUTPUT ${idSample}_${status}.md.bam
  """
}

// Creating a TSV file to restart from this step
markDuplicatesTSV.map { idPatient, status, idSample, bam, bai ->
  gender = patientGenders[idPatient]
  "${idPatient}\t${gender}\t${status}\t${idSample}\t${directoryMap.duplicateMarked}/${bam}\t${directoryMap.duplicateMarked}/${bai}\n"
}.collectFile(
  name: 'duplicateMarked.tsv', sort: true, storeDir: directoryMap.duplicateMarked
)

duplicateMarkedBams = duplicateMarkedBams.map {
    idPatient, bam, bai ->
    tag = bam.baseName.tokenize('.')[0]
    status   = tag[-1..-1].toInteger()
    idSample = tag.take(tag.length()-2)
    [idPatient, status, idSample, bam, bai]
}

if (params.verbose) duplicateMarkedBams = duplicateMarkedBams.view {
  "Realigned BAM to CreateRecalibrationTable:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  Files : [${it[3].fileName}, ${it[4].fileName}]"
}

(mdBam, mdBamToJoin) = duplicateMarkedBams.into(2)

process CreateRecalibrationTable {
  tag {idPatient + "-" + idSample}

  publishDir directoryMap.duplicateMarked, mode: params.publishDirMode, overwrite: false

  input:
    set idPatient, status, idSample, file(bam), file(bai) from mdBam // realignedBam
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
    set idPatient, status, idSample, file("${idSample}.recal.table") into recalibrationTable
    set idPatient, status, idSample, val("${idSample}_${status}.md.bam"), val("${idSample}_${status}.md.bai"), val("${idSample}.recal.table") into recalibrationTableTSV

  when: ( step == 'mapping' ) && !params.onlyQC

  script:
  known = knownIndels.collect{ "--known-sites ${it}" }.join(' ')
  """
  gatk --java-options -Xmx${task.memory.toGiga()}g \
  BaseRecalibrator \
  --input ${bam} \
  --output ${idSample}.recal.table \
  --tmp-dir /tmp \
  -R ${genomeFile} \
  -L ${intervals} \
  --known-sites ${dbsnp} \
  ${known} \
  --verbosity INFO
  """
}

// Create a TSV file to restart from this step
recalibrationTableTSV.map { idPatient, status, idSample, bam, bai, recalTable ->
  gender = patientGenders[idPatient]
  "${idPatient}\t${gender}\t${status}\t${idSample}\t${directoryMap.duplicateMarked}/${bam}\t${directoryMap.duplicateMarked}/${bai}\t${directoryMap.duplicateMarked}/${recalTable}\n"
}.collectFile(
  name: 'duplicateMarked.tsv', sort: true, storeDir: directoryMap.duplicateMarked
)

recalibrationTable = mdBamToJoin.join(recalibrationTable, by:[0,1,2])

if (step == 'recalibrate') recalibrationTable = bamFiles

if (params.verbose) recalibrationTable = recalibrationTable.view {
  "Base recalibrated table for RecalibrateBam:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  Files : [${it[3].fileName}, ${it[4].fileName}, ${it[5].fileName}]"
}

process RecalibrateBam {
  tag {idPatient + "-" + idSample}

  publishDir directoryMap.recalibrated, mode: params.publishDirMode

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

  when: !params.onlyQC

  script:
  """
  gatk --java-options -Xmx${task.memory.toGiga()}g \
  ApplyBQSR \
  -R ${genomeFile} \
  --input ${bam} \
  --output ${idSample}.recal.bam \
  -L ${intervals} \
  --create-output-bam-index true \
  --bqsr-recal-file ${recalibrationReport}
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

// Remove recalTable from Channels to match inputs for Process to avoid:
// WARN: Input tuple does not match input set cardinality declared by process...
(bamForBamQC, bamForSamToolsStats) = recalibratedBamForStats.map{ it[0..4] }.into(2)

process RunSamtoolsStats {
  tag {idPatient + "-" + idSample}

  publishDir directoryMap.samtoolsStats, mode: params.publishDirMode

  input:
    set idPatient, status, idSample, file(bam), file(bai) from bamForSamToolsStats

  output:
    file ("${bam}.samtools.stats.out") into samtoolsStatsReport

  when: !params.noReports

  script: QC.samtoolsStats(bam)
}

if (params.verbose) samtoolsStatsReport = samtoolsStatsReport.view {
  "SAMTools stats report:\n\
  File  : [${it.fileName}]"
}

process RunBamQCrecalibrated {
  tag {idPatient + "-" + idSample}

  publishDir directoryMap.bamQC, mode: params.publishDirMode

  input:
    set idPatient, status, idSample, file(bam), file(bai) from bamForBamQC

  output:
    file(idSample) into bamQCrecalibratedReport

  when: !params.noReports && !params.noBAMQC

  script:
  """
  qualimap --java-mem-size=${task.memory.toGiga()}G \
  bamqc \
  -bam ${bam} \
  --paint-chromosome-limits \
  --genome-gc-distr HUMAN \
  -nt ${task.cpus} \
  -skip-duplicated \
  --skip-dup-mode 0 \
  -outdir ${idSample} \
  -outformat HTML
  """
}

if (params.verbose) bamQCrecalibratedReport = bamQCrecalibratedReport.view {
  "BamQC report:\n\
  Dir   : [${it.fileName}]"
}

/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def checkParamReturnFile(item) {
  params."${item}" = params.genomes[params.genome]."${item}"
  return file(params."${item}")
}

def checkUppmaxProject() {
  // check if UPPMAX project number is specified
  return !(workflow.profile == 'slurm' && !params.project)
}

def checkExactlyOne(list) {
  def n = 0
  list.each{n += it ? 1 : 0}
  return n == 1
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
    'recalibrate'
  ]
}

def extractFastq(tsvFile) {
  // Channeling the TSV file containing FASTQ.
  // Format is: "subject gender status sample lane fastq1 fastq2"
  Channel.from(tsvFile)
  .splitCsv(sep: '\t')
  .map { row ->
    SarekUtils.checkNumberOfItem(row, 7)
    def idPatient  = row[0]
    def gender     = row[1]
    def status     = SarekUtils.returnStatus(row[2].toInteger())
    def idSample   = row[3]
    def idRun      = row[4]
    def fastqFile1 = SarekUtils.returnFile(row[5])
    def fastqFile2 = SarekUtils.returnFile(row[6])

    SarekUtils.checkFileExtension(fastqFile1,".fastq.gz")
    SarekUtils.checkFileExtension(fastqFile2,".fastq.gz")

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
  Channel.from(tsvFile)
    .splitCsv(sep: '\t')
    .map { row ->
    SarekUtils.checkNumberOfItem(row, 7)
    def idPatient  = row[0]
    def gender     = row[1]
    def status     = SarekUtils.returnStatus(row[2].toInteger())
    def idSample   = row[3]
    def bamFile    = SarekUtils.returnFile(row[4])
    def baiFile    = SarekUtils.returnFile(row[5])
    def recalTable = SarekUtils.returnFile(row[6])

    SarekUtils.checkFileExtension(bamFile,".bam")
    SarekUtils.checkFileExtension(baiFile,".bai")
    SarekUtils.checkFileExtension(recalTable,".recal.table")

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
  log.info "Cont Engine : " + workflow.containerEngine
  log.info "Out Dir     : " + params.outDir
  log.info "TSV file    : ${tsvFile}"
  log.info "Genome      : " + params.genome
  log.info "Genome_base : " + params.genome_base
  log.info "Step        : " + step
  log.info "Containers"
  if (params.repository != "") log.info "  Repository   : " + params.repository
  if (params.containerPath != "") log.info "  ContainerPath: " + params.containerPath
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

def sarekMessage() {
  // Display Sarek message
  log.info "Sarek - Workflow For Somatic And Germline Variations ~ ${workflow.manifest.version} - " + this.grabRevision() + (workflow.commitId ? " [${workflow.commitId}]" : "")
}

def startMessage() {
  // Display start message
  SarekUtils.sarek_ascii()
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
