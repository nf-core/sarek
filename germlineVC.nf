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
 - CreateIntervalBeds - Create and sort intervals into bed files
 - RunHaplotypecaller - Run HaplotypeCaller for Germline Variant Calling (Parallelized processes)
 - RunGenotypeGVCFs - Run HaplotypeCaller for Germline Variant Calling (Parallelized processes)
 - ConcatVCF - Merge results from paralellized callers
 - RunSingleStrelka - Run Strelka for Germline Variant Calling
 - RunSingleManta - Run Manta for Single Structural Variant Calling
 - RunBcftoolsStats - Run BCFTools stats on vcf files
 - RunVcftools - Run VCFTools on vcf files
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

if (params.help) exit 0, helpMessage()
if (!SarekUtils.isAllowedParams(params)) exit 1, "params unknown, see --help for more information"
if (!checkUppmaxProject()) exit 1, "No UPPMAX project ID found! Use --project <UPPMAX Project ID>"

// Check for awsbatch profile configuration
// make sure queue is defined
if (workflow.profile == 'awsbatch') {
    if (!params.awsqueue) exit 1, "Provide the job queue for aws batch!"
}

tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase()} : []

toolList = defineToolList()
if (!SarekUtils.checkParameterList(tools,toolList)) exit 1, 'Unknown tool(s), see --help for more information'

referenceMap = defineReferenceMap(tools)
if (!SarekUtils.checkReferenceMap(referenceMap)) exit 1, 'Missing Reference file(s), see --help for more information'

if (params.test && params.genome in ['GRCh37', 'GRCh38']) {
  referenceMap.intervals = file("$workflow.projectDir/repeats/tiny_${params.genome}.list")
}

// TODO
// FreeBayes does not need recalibrated BAMs, but we need to test whether
// the channels are set up correctly when we disable it

tsvPath = ''
if (params.sample) tsvPath = params.sample
else tsvPath = "${params.outDir}/Preprocessing/Recalibrated/recalibrated.tsv"

// Set up the bamFiles channel

bamFiles = Channel.empty()
if (tsvPath) {
  tsvFile = file(tsvPath)
  bamFiles = SarekUtils.extractBams(tsvFile, "germline")
} else exit 1, 'No sample were defined, see --help'

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

startMessage()

if (params.verbose) bamFiles = bamFiles.view {
  "BAMs to process:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  Files : [${it[3].fileName}, ${it[4].fileName}]"
}

// Here we have a recalibrated bam set
// The sample tsv config file which is formatted like: "idPatient status idSample bamFile baiFile"
// Manta will be run in Germline mode, the Tumor mode is run in somaticVC.nf
// HaplotypeCaller and Strelka will be run for Normal and Tumor samples

(bamsForGermlineManta, bamsForGermlineStrelka, bamFiles) = bamFiles.into(3)

// To speed Variant Callers up we are chopping the reference into smaller pieces.
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
        t = (\$3 - \$2) / ${params.nucleotidesPerSecond}
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
    }' ${intervals}
    """
  else
    """
    awk -vFS="[:-]" '{
      name = sprintf("%s_%d-%d", \$1, \$2, \$3);
      printf("%s\\t%d\\t%d\\n", \$1, \$2-1, \$3) > name ".bed"
    }' ${intervals}
    """
}

bedIntervals = bedIntervals
  .map { intervalFile ->
    def duration = 0.0
    for (line in intervalFile.readLines()) {
      final fields = line.split('\t')
      if (fields.size() >= 5) duration += fields[4].toFloat()
      else {
        start = fields[1].toInteger()
        end = fields[2].toInteger()
        duration += (end - start) / params.nucleotidesPerSecond
      }
    }
    [duration, intervalFile]
  }.toSortedList({ a, b -> b[0] <=> a[0] })
  .flatten().collate(2)
  .map{duration, intervalFile -> intervalFile}

if (params.verbose) bedIntervals = bedIntervals.view {
  "  Interv: ${it.baseName}"
}

bamsForHC = bamFiles.combine(bedIntervals)

process RunHaplotypecaller {
  tag {idSample + "-" + intervalBed.baseName}

  input:
    set idPatient, status, idSample, file(bam), file(bai), file(intervalBed) from bamsForHC
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(dbsnp), file(dbsnpIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.dbsnp,
      referenceMap.dbsnpIndex
    ])

  output:
    set val("HaplotypeCallerGVCF"), idPatient, idSample, file("${intervalBed.baseName}_${idSample}.g.vcf") into hcGenomicVCF
    set idPatient, status, idSample, file(intervalBed), file("${intervalBed.baseName}_${idSample}.g.vcf") into vcfsToGenotype

  when: 'haplotypecaller' in tools && !params.onlyQC

  script:
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
    HaplotypeCaller \
    -R ${genomeFile} \
    -I ${bam} \
    -L ${intervalBed} \
    -D ${dbsnp} \
    -O ${intervalBed.baseName}_${idSample}.g.vcf \
    -ERC GVCF
  """
}
hcGenomicVCF = hcGenomicVCF.groupTuple(by:[0,1,2])

if (params.noGVCF) hcGenomicVCF.close()

process RunGenotypeGVCFs {
  tag {idSample + "-" + intervalBed.baseName}

  input:
    set idPatient, status, idSample, file(intervalBed), file(gvcf) from vcfsToGenotype
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(dbsnp), file(dbsnpIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.dbsnp,
      referenceMap.dbsnpIndex
    ])

  output:
    set val("HaplotypeCaller"), idPatient, idSample, file("${intervalBed.baseName}_${idSample}.vcf") into hcGenotypedVCF

  when: 'haplotypecaller' in tools && !params.onlyQC

  script:
  // Using -L is important for speed and we have to index the interval files also
  """
  gatk --java-options -Xmx${task.memory.toGiga()}g \
    IndexFeatureFile -F ${gvcf}

  gatk --java-options -Xmx${task.memory.toGiga()}g \
    GenotypeGVCFs \
    -R ${genomeFile} \
    -L ${intervalBed} \
    -D ${dbsnp} \
    -V ${gvcf} \
    -O ${intervalBed.baseName}_${idSample}.vcf
  """
}
hcGenotypedVCF = hcGenotypedVCF.groupTuple(by:[0,1,2,3])

// we are merging the VCFs that are called separatelly for different intervals
// so we can have a single sorted VCF containing all the calls for a given caller

vcfsToMerge = hcGenomicVCF.mix(hcGenotypedVCF)
if (params.verbose) vcfsToMerge = vcfsToMerge.view {
  "VCFs To be merged:\n\
  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: [${it[2]}]\n\
  Files : ${it[3].fileName}"
}

process ConcatVCF {
  tag {variantCaller + "-" + idSample}

  publishDir "${params.outDir}/VariantCalling/${idPatient}/${"$variantCaller"}", mode: params.publishDirMode

  input:
    set variantCaller, idPatient, idSample, file(vcFiles) from vcfsToMerge
    file(genomeIndex) from Channel.value(referenceMap.genomeIndex)
    file(targetBED) from Channel.value(params.targetBED ? file(params.targetBED) : "null")

  output:
    // we have this funny *_* pattern to avoid copying the raw calls to publishdir
    set variantCaller, idPatient, idSample, file("*_*.vcf.gz"), file("*_*.vcf.gz.tbi") into vcfConcatenated

  when: 'haplotypecaller' in tools && !params.onlyQC

  script:
  if (variantCaller == 'HaplotypeCaller') outputFile = "${variantCaller}_${idSample}.vcf"
  else if (variantCaller == 'HaplotypeCallerGVCF') outputFile = "haplotypecaller_${idSample}.g.vcf"
  options = params.targetBED ? "-t ${targetBED}" : ""
  """
  concatenateVCFs.sh -i ${genomeIndex} -c ${task.cpus} -o ${outputFile} ${options}
  """
}

if (params.verbose) vcfConcatenated = vcfConcatenated.view {
  "Variant Calling output:\n\
  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: ${it[2]}\n\
  Files : ${it[3].fileName}\n\
  Index : ${it[4].fileName}"
}

process RunSingleStrelka {
  tag {idSample}

  publishDir "${params.outDir}/VariantCalling/${idPatient}/Strelka", mode: params.publishDirMode

  input:
    set idPatient, status, idSample, file(bam), file(bai) from bamsForGermlineStrelka
    file(targetBED) from Channel.value(params.targetBED ? file(params.targetBED) : "null")
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex
    ])

  output:
    set val("Strelka"), idPatient, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into singleStrelkaOutput

  when: 'strelka' in tools && !params.onlyQC

  script:
  beforeScript = params.targetBED ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
  options = params.targetBED ? "--exome --callRegions call_targets.bed.gz" : ""
  """
  ${beforeScript}
  configureStrelkaGermlineWorkflow.py \
  --bam ${bam} \
  --referenceFasta ${genomeFile} \
  ${options} \
  --runDir Strelka

  python Strelka/runWorkflow.py -m local -j ${task.cpus}
  mv Strelka/results/variants/genome.*.vcf.gz Strelka_${idSample}_genome.vcf.gz
  mv Strelka/results/variants/genome.*.vcf.gz.tbi Strelka_${idSample}_genome.vcf.gz.tbi
  mv Strelka/results/variants/variants.vcf.gz Strelka_${idSample}_variants.vcf.gz
  mv Strelka/results/variants/variants.vcf.gz.tbi Strelka_${idSample}_variants.vcf.gz.tbi
  """
}

if (params.verbose) singleStrelkaOutput = singleStrelkaOutput.view {
  "Variant Calling output:\n\
  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: ${it[2]}\n\
  Files : ${it[3].fileName}\n\
  Index : ${it[4].fileName}"
}

process RunSingleManta {
  tag {idSample + " - Single Diploid"}

  publishDir "${params.outDir}/VariantCalling/${idPatient}/Manta", mode: params.publishDirMode

  input:
    set idPatient, status, idSample, file(bam), file(bai) from bamsForGermlineManta
    file(targetBED) from Channel.value(params.targetBED ? file(params.targetBED) : "null")
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex
    ])

  output:
    set val("Manta"), idPatient, idSample,  file("*.vcf.gz"), file("*.vcf.gz.tbi") into singleMantaOutput

  when: 'manta' in tools  && status == 0 && !params.onlyQC

  script:
  beforeScript = params.targetBED ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
  options = params.targetBED ? "--exome --callRegions call_targets.bed.gz" : ""
  """
  ${beforeScript}
  configManta.py \
  --bam ${bam} \
  --reference ${genomeFile} \
  ${options} \
  --runDir Manta

  python Manta/runWorkflow.py -m local -j ${task.cpus}

  mv Manta/results/variants/candidateSmallIndels.vcf.gz \
    Manta_${idSample}.candidateSmallIndels.vcf.gz
  mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
    Manta_${idSample}.candidateSmallIndels.vcf.gz.tbi
  mv Manta/results/variants/candidateSV.vcf.gz \
    Manta_${idSample}.candidateSV.vcf.gz
  mv Manta/results/variants/candidateSV.vcf.gz.tbi \
    Manta_${idSample}.candidateSV.vcf.gz.tbi
  mv Manta/results/variants/diploidSV.vcf.gz \
    Manta_${idSample}.diploidSV.vcf.gz
  mv Manta/results/variants/diploidSV.vcf.gz.tbi \
    Manta_${idSample}.diploidSV.vcf.gz.tbi
  """
}

if (params.verbose) singleMantaOutput = singleMantaOutput.view {
  "Variant Calling output:\n\
  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: ${it[2]}\n\
  Files : ${it[3].fileName}\n\
  Index : ${it[4].fileName}"
}

vcfForQC = Channel.empty().mix(
  vcfConcatenated.map {
    variantcaller, idPatient, idSample, vcf, tbi ->
    [variantcaller, vcf]
  },
  singleStrelkaOutput.map {
    variantcaller, idPatient, idSample, vcf, tbi ->
    [variantcaller, vcf[1]]
  },
  singleMantaOutput.map {
    variantcaller, idPatient, idSample, vcf, tbi ->
    [variantcaller, vcf[2]]
  })

(vcfForBCFtools, vcfForVCFtools) = vcfForQC.into(2)

process RunBcftoolsStats {
  tag {"${variantCaller} - ${vcf}"}

  publishDir "${params.outDir}/Reports/BCFToolsStats", mode: params.publishDirMode

  input:
    set variantCaller, file(vcf) from vcfForBCFtools

  output:
    file ("*.bcf.tools.stats.out") into bcfReport

  when: !params.noReports

  script: QC.bcftools(vcf)
}

if (params.verbose) bcfReport = bcfReport.view {
  "BCFTools stats report:\n\
  File  : [${it.fileName}]"
}

bcfReport.close()

process RunVcftools {
  tag {"${variantCaller} - ${vcf}"}

  publishDir "${params.outDir}/Reports/VCFTools", mode: params.publishDirMode

  input:
    set variantCaller, file(vcf) from vcfForVCFtools

  output:
    file ("${reducedVCF}.*") into vcfReport

  when: !params.noReports

  script:
    reducedVCF = SarekUtils.reduceVCF(vcf)
    QC.vcftools(vcf)
}

if (params.verbose) vcfReport = vcfReport.view {
  "VCFTools stats report:\n\
  File  : [${it.fileName}]"
}

vcfReport.close()

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

def defineReferenceMap(tools) {
  if (!(params.genome in params.genomes)) exit 1, "Genome ${params.genome} not found in configuration"
  def referenceMap =
  [
    'genomeFile'       : checkParamReturnFile("genomeFile"),
    'genomeIndex'      : checkParamReturnFile("genomeIndex"),
    'intervals'        : checkParamReturnFile("intervals")
  ]
  if ('haplotypecaller' in tools) {
    referenceMap.putAll(
      'dbsnp'            : checkParamReturnFile("dbsnp"),
      'dbsnpIndex'       : checkParamReturnFile("dbsnpIndex"),
      'genomeDict'       : checkParamReturnFile("genomeDict")
    )
  }
  return referenceMap
}

def defineToolList() {
  return [
    'ascat',
    'freebayes',
    'haplotypecaller',
    'manta',
    'mutect2',
    'strelka'
  ]
}

def grabRevision() {
  // Return the same string executed from github or not
  return workflow.revision ?: workflow.commitId ?: workflow.scriptId.substring(0,10)
}

def helpMessage() {
  // Display help message
  this.sarekMessage()
  log.info "    Usage:"
  log.info "       nextflow run germlineVC.nf --sample <file.tsv> [--tools TOOL[,TOOL]] --genome <Genome>"
  log.info "    --sample <file.tsv>"
  log.info "       Specify a TSV file containing paths to sample files."
  log.info "    --test"
  log.info "       Use a test sample."
  log.info "    --noReports"
  log.info "       Disable QC tools and MultiQC to generate a HTML report"
  log.info "    --tools"
  log.info "       Option to configure which tools to use in the workflow."
  log.info "         Different tools to be separated by commas."
  log.info "       Possible values are:"
  log.info "         strelka (use Strelka for VC)"
  log.info "         haplotypecaller (use HaplotypeCaller for normal bams VC)"
  log.info "         manta (use Manta for SV)"
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
  log.info "Out Dir     : " + params.outDir
  log.info "TSV file    : " + tsvFile
  log.info "Genome      : " + params.genome
  log.info "Genome_base : " + params.genome_base
  log.info "Target BED  : " + params.targetBED
  log.info "Tools       : " + tools.join(', ')
  log.info "Containers"
  if (params.repository != "") log.info "  Repository   : " + params.repository
  if (params.containerPath != "") log.info "  ContainerPath: " + params.containerPath
  log.info "Reference files used:"
  log.info "  dbsnp       :\n\t" + referenceMap.dbsnp
  log.info "\t" + referenceMap.dbsnpIndex
  log.info "  genome      :\n\t" + referenceMap.genomeFile
  log.info "\t" + referenceMap.genomeDict
  log.info "\t" + referenceMap.genomeIndex
  log.info "  intervals   :\n\t" + referenceMap.intervals
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
