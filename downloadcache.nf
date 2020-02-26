#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/sarek
========================================================================================
 nf-core/sarek Analysis Pipeline.
 @Homepage
 https://sarek.scilifelab.se/
 @Documentation
 https://github.com/nf-core/sarek/README.md
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

Usage:
    --help
      you're reading it

DOWNLOAD CACHE:
  nextflow run build.nf --download_cache [--snpEff_cache <pathToSNPEFFcache>] [--vep_cache <pathToVEPcache>]
                                         [--cadd_cache <pathToCADDcache> --cadd_version <CADD Version>]
    --download_cache
      Will download specified cache
    --snpEff_cache <Directoy>
      Specify path to snpEff cache
      If none, will use snpEff version specified in configuration
      Will use snpEff cache version for ${params.genome}: ${params.genomes[params.genome].snpeffDb} in igenomes configuration file:
      Change with --genome or in configuration files
    --vep_cache <Directoy>
      Specify path to VEP cache
      If none, will use VEP version specified in configuration
      Will use VEP cache version for ${params.genome}: ${params.genomes[params.genome].vepCacheVersion} in igenomes configuration file:
      Change with --genome or in configuration files
    --cadd_cache <Directoy>
      Specify path to CADD cache
      Will use CADD version specified
    --cadd_version <version>
      Will specify which CADD version to download
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help) exit 0, helpMessage()

// Default value for params
params.build = null
params.offline = null
params.cadd_cache = null
params.cadd_version = 'v1.5'
params.genome = 'GRCh37'
params.snpEff_cache = null
params.vep_cache = null

ch_referencesFiles = Channel.empty()

pathToSource = params.offline ? "data/reference/" : "https://github.com/nf-core/test-datasets/raw/sarek/reference"

if (params.build) ch_referencesFiles = ch_referencesFiles.mix(
  Channel.fromPath("${pathToSource}/1000G_phase1.indels.b37.small.vcf.gz"),
  Channel.fromPath("${pathToSource}/1000G_phase3_20130502_SNP_maf0.3.small.loci"),
  Channel.fromPath("${pathToSource}/1000G_phase3_20130502_SNP_maf0.3.small.loci.gc"),
  Channel.fromPath("${pathToSource}/Mills_and_1000G_gold_standard.indels.b37.small.vcf.gz"),
  Channel.fromPath("${pathToSource}/dbsnp_138.b37.small.vcf.gz"),
  Channel.fromPath("${pathToSource}/human_g1k_v37_decoy.small.fasta.gz"),
  Channel.fromPath("${pathToSource}/small.intervals"))

ch_referencesFiles = ch_referencesFiles.dump(tag:'Reference Files')

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if ( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

if ( workflow.profile == 'awsbatch') {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  // Check outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
  // Prevent trace files to be stored on S3 since S3 does not support rolling files.
  if (workflow.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile == 'awsbatch'){
   summary['AWS Region']    = params.awsregion
   summary['AWS Queue']     = params.awsqueue
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if (params.email) {
  summary['E-mail Address']  = params.email
  summary['MultiQC maxsize'] = params.maxMultiqcEmailFileSize
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
if (params.monochrome_logs) log.info "----------------------------------------------------"
else log.info "\033[2m----------------------------------------------------\033[0m"

// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-sarek-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/sarek Workflow Summary'
    section_href: 'https://github.com/nf-core/sarek'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}

/*
================================================================================
=                   D  O  W  N  L  O  A  D      C  A  C  H  E                  =
================================================================================
*/

process BuildCache_snpEff {
  tag {snpeffDb}

  publishDir params.snpEff_cache, mode: params.publishDirMode

  input:
    val snpeffDb from Channel.value(params.genomes[params.genome].snpeffDb)

  output:
    file("*")

  when: params.snpEff_cache && params.download_cache && !params.offline

  script:
  """
  snpEff download -v ${snpeffDb} -dataDir \${PWD}
  """
}

process BuildCache_VEP {
  tag {"${species}_${cache_version}_${genome}"}

  publishDir "${params.vep_cache}/${species}", mode: params.publishDirMode

  input:
    val cache_version from Channel.value(params.genomes[params.genome].vepCacheVersion)

  output:
    file("*")

  when: params.vep_cache && params.download_cache && !params.offline

  script:
  genome = params.genome
  species = genome =~ "GRCh3*" ? "homo_sapiens" : genome =~ "GRCm3*" ? "mus_musculus" : ""
  """
  vep_install \
    -a cf \
    -c . \
    -s ${species} \
    -v ${cache_version} \
    -y ${genome} \
    --CACHE_VERSION ${cache_version} \
    --CONVERT \
    --NO_HTSLIB --NO_TEST --NO_BIOPERL --NO_UPDATE

  mv ${species}/* .
  rm -rf ${species}
  """
}

caddFileToDownload = (params.cadd_version) && (params.genome == "GRCh37" || params.genome == "GRCh38") ?
  Channel.from(
    "https://krishna.gs.washington.edu/download/CADD/${params.cadd_version}/${params.genome}/InDels_inclAnno.tsv.gz",
    "https://krishna.gs.washington.edu/download/CADD/${params.cadd_version}/${params.genome}/whole_genome_SNVs_inclAnno.tsv.gz"
  ) : Channel.empty()

process DownloadCADD {
  tag {caddFile}

  publishDir "${params.cadd_cache}/${params.genome}", mode: params.publishDirMode

  input:
    val(caddFile) from caddFileToDownload

  output:
    set file("*.tsv.gz"), file("*.tsv.gz.tbi")

  when: params.cadd_cache && params.download_cache && !params.offline

  script:
  """
  wget --quiet ${caddFile}
  wget --quiet ${caddFile}.tbi
  """
}

def nfcoreHeader(){
    // Log colors ANSI codes
    c_reset  = params.monochrome_logs ? '' : "\033[0m";
    c_dim    = params.monochrome_logs ? '' : "\033[2m";
    c_black  = params.monochrome_logs ? '' : "\033[0;30m";
    c_red    = params.monochrome_logs ? '' : "\033[0;31m";
    c_green  = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue   = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan   = params.monochrome_logs ? '' : "\033[0;36m";
    c_white  = params.monochrome_logs ? '' : "\033[0;37m";

    return """    ${c_dim}----------------------------------------------------${c_reset}
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
        ${c_white}____${c_reset}
      ${c_white}.´ _  `.${c_reset}
     ${c_white}/  ${c_green}|\\${c_reset}`-_ \\${c_reset}     ${c_blue} __        __   ___     ${c_reset}
    ${c_white}|   ${c_green}| \\${c_reset}  `-|${c_reset}    ${c_blue}|__`  /\\  |__) |__  |__/${c_reset}
     ${c_white}\\ ${c_green}|   \\${c_reset}  /${c_reset}     ${c_blue}.__| /¯¯\\ |  \\ |___ |  \\${c_reset}
      ${c_white}`${c_green}|${c_reset}____${c_green}\\${c_reset}´${c_reset}

    ${c_purple}  nf-core/sarek v${workflow.manifest.version}${c_reset}
    ${c_dim}----------------------------------------------------${c_reset}
    """.stripIndent()
}

def checkHostname(){
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames){
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)){
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}

/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def checkFile(it) {
  // Check file existence
  final f = file(it)
  if (!f.exists()) exit 1, "Missing file: ${it}, see --help for more information"
  return true
}
