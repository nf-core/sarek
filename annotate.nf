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
 - RunBcftoolsStats - Run BCFTools stats on vcf files
 - RunVcftools - Run VCFTools on vcf files
 - RunSnpeff - Run snpEff for annotation of vcf files
 - RunVEP - Run VEP for annotation of vcf files
 - CompressVCF - Compress and index vcf files using tabix
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

annotateTools = params.annotateTools ? params.annotateTools.split(',').collect{it.trim().toLowerCase()} : []
annotateVCF = params.annotateVCF ? params.annotateVCF.split(',').collect{it.trim()} : []
tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase()} : []

toolList = defineToolList()

if (!SarekUtils.checkParameterList(tools,toolList)) exit 1, 'Unknown tool(s), see --help for more information'

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

startMessage()

vcfToAnnotate = Channel.create()
vcfNotToAnnotate = Channel.create()

if (annotateVCF == []) {
// Sarek, by default, annotates all available vcfs that it can find in the VariantCalling directory
// Excluding vcfs from FreeBayes, and g.vcf from HaplotypeCaller
// Basically it's: VariantCalling/*/{HaplotypeCaller,Manta,MuTect2,Strelka}/*.vcf.gz
// Without *SmallIndels.vcf.gz from Manta, and *.genome.vcf.gz from Strelka
// The small snippet `vcf.minus(vcf.fileName)[-2]` catches idPatient
// This field is used to output final annotated VCFs in the correct directory
  Channel.empty().mix(
    Channel.fromPath("${params.outDir}/VariantCalling/*/HaplotypeCaller/*.vcf.gz")
      .flatten().map{vcf -> ['haplotypecaller', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
    Channel.fromPath("${params.outDir}/VariantCalling/*/Manta/*[!candidate]SV.vcf.gz")
      .flatten().map{vcf -> ['manta', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
    Channel.fromPath("${params.outDir}/VariantCalling/*/MuTect2/*.vcf.gz")
      .flatten().map{vcf -> ['mutect2', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
    Channel.fromPath("${params.outDir}/VariantCalling/*/Strelka/*{somatic,variant}*.vcf.gz")
      .flatten().map{vcf -> ['strelka', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
  ).choice(vcfToAnnotate, vcfNotToAnnotate) {
    annotateTools == [] || (annotateTools != [] && it[0] in annotateTools) ? 0 : 1
  }
} else if (annotateTools == []) {
// Annotate user-submitted VCFs
// If user-submitted, Sarek assume that the idPatient should be assumed automatically
  vcfToAnnotate = Channel.fromPath(annotateVCF)
    .map{vcf -> ['userspecified', vcf.minus(vcf.fileName)[-2].toString(), vcf]}
} else exit 1, "specify only tools or files to annotate, not both"

vcfNotToAnnotate.close()

// as now have the list of VCFs to annotate, the first step is to annotate with allele frequencies, if there are any

(vcfForBCFtools, vcfForVCFtools, vcfForSnpeff, vcfForVep) = vcfToAnnotate.into(4)

vcfForVep = vcfForVep.map {
  variantCaller, idPatient, vcf ->
  ["VEP", variantCaller, idPatient, vcf, null]
}

process RunBcftoolsStats {
  tag {"${idPatient} - ${variantCaller} - ${vcf}"}

  publishDir "${params.outDir}/Reports/BCFToolsStats", mode: params.publishDirMode

  input:
    set variantCaller, idPatient, file(vcf) from vcfForBCFtools

  output:
    file ("*.bcf.tools.stats.out") into bcfReport

  when: !params.noReports

  script: QC.bcftools(vcf)
}

if (params.verbose) bcfReport = bcfReport.view {
  "BCFTools stats report:\n" +
  "File  : [${it.fileName}]"
}

process RunVcftools {
  tag {"${idPatient} - ${variantCaller} - ${vcf}"}

  publishDir "${params.outDir}/Reports/VCFTools", mode: params.publishDirMode

  input:
    set variantCaller, idPatient, file(vcf) from vcfForVCFtools

  output:
    file ("${reducedVCF}.*") into vcfReport

  when: !params.noReports

  script:
    reducedVCF = SarekUtils.reduceVCF(vcf)
    QC.vcftools(vcf)
}

if (params.verbose) vcfReport = vcfReport.view {
  "VCFTools stats report:\n" +
  "Files : [${it.fileName}]"
}

process RunSnpeff {
  tag {"${idPatient} - ${variantCaller} - ${vcf}"}

  publishDir params.outDir, mode: params.publishDirMode, saveAs: {
    if (it == "${reducedVCF}_snpEff.ann.vcf") null
    else "Annotation/${idPatient}/snpEff/${it}"
  }

  input:
    set variantCaller, idPatient, file(vcf) from vcfForSnpeff
    file dataDir from Channel.value(params.snpEff_cache ? file(params.snpEff_cache) : "null")
    val snpeffDb from Channel.value(params.genomes[params.genome].snpeffDb)

  output:
    set file("${reducedVCF}_snpEff.genes.txt"), file("${reducedVCF}_snpEff.csv"), file("${reducedVCF}_snpEff.summary.html") into snpeffOutput
    set val("snpEff"), variantCaller, idPatient, file("${reducedVCF}_snpEff.ann.vcf") into snpeffVCF

  when: 'snpeff' in tools || 'merge' in tools

  script:
  reducedVCF = SarekUtils.reduceVCF(vcf)
  cache = (params.snpEff_cache && params.annotation_cache) ? "-dataDir \${PWD}/${dataDir}" : ""
  """
  snpEff -Xmx${task.memory.toGiga()}g \
  ${snpeffDb} \
  -csvStats ${reducedVCF}_snpEff.csv \
  -nodownload \
  ${cache} \
  -canon \
  -v \
  ${vcf} \
  > ${reducedVCF}_snpEff.ann.vcf

  mv snpEff_summary.html ${reducedVCF}_snpEff.summary.html
  """
}

if (params.verbose) snpeffOutput = snpeffOutput.view {
  "snpEff report:\n" +
  "File  : ${it.fileName}"
}

if ('merge' in tools) {
  // When running in the 'merge' mode
  // snpEff output is used as VEP input
  // Used a feedback loop from vcfCompressed
  // https://github.com/nextflow-io/patterns/tree/master/feedback-loop

  vcfCompressed = Channel.create()

  vcfForVep = Channel.empty().mix(
    vcfCompressed.until({ it[0]=="merge" })
  )
}

process RunVEP {
  tag {"${idPatient} - ${variantCaller} - ${vcf}"}

  publishDir params.outDir, mode: params.publishDirMode, saveAs: {
    if (it == "${reducedVCF}_VEP.summary.html") "Annotation/${idPatient}/VEP/${it}"
    else null
  }

  input:
    set annotator, variantCaller,  idPatient, file(vcf), file(idx) from vcfForVep
    file dataDir from Channel.value(params.vep_cache ? file(params.vep_cache) : "null")
    val cache_version from Channel.value(params.genomes[params.genome].vepCacheVersion)
    set file(cadd_WG_SNVs), file(cadd_WG_SNVs_tbi), file(cadd_InDels), file(cadd_InDels_tbi) from Channel.value([
      params.cadd_WG_SNVs ? file(params.cadd_WG_SNVs) : "null",
      params.cadd_WG_SNVs_tbi ? file(params.cadd_WG_SNVs_tbi) : "null",
      params.cadd_InDels ? file(params.cadd_InDels) : "null",
      params.cadd_InDels_tbi ? file(params.cadd_InDels_tbi) : "null"
    ])

  output:
    set finalAnnotator, variantCaller, idPatient, file("${reducedVCF}_VEP.ann.vcf") into vepVCF
    file("${reducedVCF}_VEP.summary.html") into vepReport

  when: 'vep' in tools || 'merge' in tools

  script:
  reducedVCF = SarekUtils.reduceVCF(vcf)
  finalAnnotator = annotator == "snpEff" ? 'merge' : 'VEP'
  genome = params.genome == 'smallGRCh37' ? 'GRCh37' : params.genome
  dir_cache = (params.vep_cache && params.annotation_cache) ? " \${PWD}/${dataDir}" : "/.vep"
  cadd = (params.cadd_cache && params.cadd_WG_SNVs && params.cadd_InDels) ? "--plugin CADD,whole_genome_SNVs.tsv.gz,InDels.tsv.gz" : ""
  genesplicer = params.genesplicer ? "--plugin GeneSplicer,/opt/conda/envs/sarek-2.3/bin/genesplicer,/opt/conda/envs/sarek-2.3/share/genesplicer-1.0-1/human,context=200,tmpdir=\$PWD/${reducedVCF}" : "--offline"
  """
  mkdir ${reducedVCF}

  vep \
  -i ${vcf} \
  -o ${reducedVCF}_VEP.ann.vcf \
  --assembly ${genome} \
  ${cadd} \
  ${genesplicer} \
  --cache \
  --cache_version ${cache_version} \
  --dir_cache ${dir_cache} \
  --everything \
  --filter_common \
  --fork ${task.cpus} \
  --format vcf \
  --per_gene \
  --stats_file ${reducedVCF}_VEP.summary.html \
  --total_length \
  --vcf

  rm -rf ${reducedVCF}
  """
}

if (params.verbose) vepReport = vepReport.view {
  "VEP report:\n" +
  "Files : ${it.fileName}"
}

vcfToCompress = snpeffVCF.mix(vepVCF)

process CompressVCF {
  tag {"${idPatient} - ${annotator} - ${vcf}"}

  publishDir "${params.outDir}/Annotation/${idPatient}/${finalAnnotator}", mode: params.publishDirMode

  input:
    set annotator, variantCaller, idPatient, file(vcf) from vcfToCompress

  output:
    set annotator, variantCaller, idPatient, file("*.vcf.gz"), file("*.vcf.gz.tbi") into (vcfCompressed, vcfCompressedoutput)

  script:
  reducedVCF = SarekUtils.reduceVCF(vcf)
  finalAnnotator = annotator == "merge" ? "VEP" : annotator
  """
  bgzip < ${vcf} > ${vcf}.gz
  tabix ${vcf}.gz
  """
}

if (params.verbose) vcfCompressedoutput = vcfCompressedoutput.view {
  "${it[2]} - ${it[0]} VCF:\n" +
  "File  : ${it[3].fileName}\n" +
  "Index : ${it[4].fileName}"
}

/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def checkUppmaxProject() {
  // check if UPPMAX project number is specified
  return !(workflow.profile == 'slurm' && !params.project)
}

def defineToolList() {
  return [
    'merge',
    'snpeff',
    'vep'
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
  log.info "       nextflow run annotate.nf --test [--step STEP] [--tools TOOL[,TOOL]] --genome <Genome>"
  log.info "    --noReports"
  log.info "       Disable QC tools and MultiQC to generate a HTML report"
  log.info "    --tools"
  log.info "       Option to configure which tools to use in the workflow."
  log.info "         Different tools to be separated by commas."
  log.info "       Possible values are:"
  log.info "         snpeff (use snpEff for Annotation of Variants)"
  log.info "         vep (use VEP for Annotation of Variants)"
  log.info "         merge (first snpEff, then feed its output VCFs to VEP)"
  log.info "    --annotateTools"
  log.info "       Option to configure which tools to annotate."
  log.info "         Different tools to be separated by commas."
  log.info "       Possible values are:"
  log.info "         haplotypecaller (Annotate HaplotypeCaller output)"
  log.info "         manta (Annotate Manta output)"
  log.info "         mutect2 (Annotate MuTect2 output)"
  log.info "         strelka (Annotate Strelka output)"
  log.info "    --annotateVCF"
  log.info "       Option to configure which vcf to annotate."
  log.info "         Different vcf to be separated by commas."
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
  log.info "Genome      : " + params.genome
  log.info "Genome_base : " + params.genome_base
  if (tools) log.info "Tools       : " + tools.join(', ')
  if (annotateTools) log.info "Annotate on : " + annotateTools.join(', ')
  if (annotateVCF) log.info "VCF files   : " +annotateVCF.join(',\n    ')
  log.info "Containers"
  if (params.repository != "") log.info "  Repository   : " + params.repository
  if (params.containerPath != "") log.info "  ContainerPath: " + params.containerPath
  log.info "  Tag          : " + params.tag
  log.info "Reference files used:"
  log.info "  snpEff DB    :\n\t" + params.genomes[params.genome].snpeffDb
  log.info "  VEP Cache    :\n\t" + params.genomes[params.genome].vepCacheVersion
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
