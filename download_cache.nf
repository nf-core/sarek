#!/usr/bin/env nextflow

/*
================================================================================
                                  nf-core/sarek
================================================================================
Started March 2016.
Ported to nf-core May 2019.
--------------------------------------------------------------------------------
nf-core/sarek:
  An open-source analysis pipeline to detect germline or somatic variants
  from whole genome or targeted sequencing
--------------------------------------------------------------------------------
 @Homepage
 https://nf-co.re/sarek
--------------------------------------------------------------------------------
 @Documentation
 https://nf-co.re/sarek/docs
--------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/sarek/download_cache.nf -profile docker --genome <genome> --help
                                      [--snpeff_cache <pathToSNPEFFcache> --snpeff_db_version <snpEff DB version>]
                                      [--vep_cache <pathToVEPcache> --vep_cache_version <VEP cache version> --species <species>]
                                      [--cadd_cache <pathToCADDcache> --cadd_version <CADD Version>]

    Options:
      --help                   [bool] You're reading it
      --snpeff_cache           [file] Path to snpEff cache
      --snpeff_db_version       [str] snpEff DB version
                                      Default: ${params.genomes[params.genome].snpeff_db}
      --vep_cache              [file] Path to VEP cache
      --vep_cache_version       [int] VEP cache version
                                      Default: ${params.genomes[params.genome].vep_cache_version}
      --species                 [str] Species
                                      Default: ${params.genomes[params.genome].species}
      --cadd_cache             [file] Path to CADD cache
      --cadd_version            [str] CADD version to download
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help) exit 0, helpMessage()

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

params.snpeff_db = params.genome ? params.genomes[params.genome].snpeff_db ?: null : null
params.species = params.genome ? params.genomes[params.genome].species ?: null : null
params.vep_cache_version = params.genome ? params.genomes[params.genome].vep_cache_version ?: null : null

ch_snpeff_db = params.snpeff_db ? Channel.value(params.snpeff_db) : "null"
ch_vep_cache_version = params.vep_cache_version ? Channel.value(params.vep_cache_version) : "null"

// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
if (params.snpeff_db)               summary['snpeffDb']              = params.snpeff_db
if (params.vep_cache_version)       summary['vepCacheVersion']       = params.vep_cache_version
if (params.species)                 summary['species']               = params.species
if (params.snpeff_cache)            summary['snpEff_cache']          = params.snpeff_cache
if (params.vep_cache)               summary['vep_cache']             = params.vep_cache
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
  tag {snpeff_db}

  publishDir params.snpeff_cache, mode: params.publish_dir_mode

  input:
    val snpeff_db from ch_snpeff_db

  output:
    file("*") into snpeff_cache_out

  when: params.snpeff_cache

  script:
  """
  snpEff download -v ${snpeff_db} -dataDir \${PWD}
  """
}

snpeff_cache_out = snpeff_cache_out.dump(tag: 'snpeff_cache_out')

process BuildCache_VEP {
  tag {"${species}_${vep_cache_version}_${genome}"}

  publishDir "${params.vep_cache}/${species}", mode: params.publish_dir_mode

  input:
    val vep_cache_version from ch_vep_cache_version
    val species from Channel.value(params.species)

  output:
    file("*") into vep_cache_out

  when: params.vep_cache

  script:
  genome = params.genome
  """
  vep_install \
    -a cf \
    -c . \
    -s ${species} \
    -y ${genome} \
    --CACHE_VERSION ${vep_cache_version} \
    --CONVERT \
    --NO_HTSLIB --NO_TEST --NO_BIOPERL --NO_UPDATE

  mv ${species}/* .
  rm -rf ${species}
  """
}

vep_cache_out = vep_cache_out.dump(tag: 'vep_cache_out')

caddFileToDownload = (params.cadd_cache && params.cadd_version) && (params.genome == "GRCh37" || params.genome == "GRCh38") ?
  Channel.from(
    "https://krishna.gs.washington.edu/download/CADD/${params.cadd_version}/${params.genome}/InDels_inclAnno.tsv.gz",
    "https://krishna.gs.washington.edu/download/CADD/${params.cadd_version}/${params.genome}/whole_genome_SNVs_inclAnno.tsv.gz"
  ) : Channel.empty()

process DownloadCADD {
  tag {caddFile}

  publishDir "${params.cadd_cache}/${params.genome}", mode: params.publish_dir_mode

  input:
    val caddFile from caddFileToDownload

  output:
    set file("*.tsv.gz"), file("*.tsv.gz.tbi") into cadd_files

  when: (params.cadd_cache && params.cadd_version) && (params.genome == "GRCh37" || params.genome == "GRCh38")

  script:
  """
  wget --quiet ${caddFile}
  wget --quiet ${caddFile}.tbi
  """
}

cadd_files = cadd_files.dump(tag: 'cadd_files')

def nfcoreHeader() {
    // Log colors ANSI codes
    c_black  = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue   = params.monochrome_logs ? '' : "\033[0;34m";
    c_dim    = params.monochrome_logs ? '' : "\033[2m";
    c_green  = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset  = params.monochrome_logs ? '' : "\033[0m";
    c_white  = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
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
    -${c_dim}--------------------------------------------------${c_reset}-
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