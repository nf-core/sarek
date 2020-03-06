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
    --help
      you're reading it

  nextflow run create_pon.nf --input sample.tsv -profile docker

    -profile                  [str] Configuration profile to use
                                    Can use multiple (comma separated)
                                    Available: conda, docker, singularity, test and more
    --input                  [file] Path to input TSV file 
    --genome                  [str] Name of iGenomes reference
    --dict                   [file] dict from the fasta reference
                                    If none provided, will be generated automatically from the fasta reference
    --fasta                  [file] fasta reference
    --fasta_fai              [file] reference index
                                    If none provided, will be generated automatically from the fasta reference
    --germline_resource      [file] Germline Resource File (AF only gnomad)
    --germline_resource_index       Germline Resource Index
                             [file] if none provided, will be generated automatically if a germlineResource file is provided
    --intervals              [file] intervals
                                    If none provided, will be generated automatically from the fasta reference
                                    Use --no_intervals to disable automatic generation
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
                         SET UP CONFIGURATION VARIABLES
================================================================================
*/

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

tsvPath = null
if (params.input && (hasExtension(params.input, "tsv") = params.input


//TODO: Define input channels following sarek logic

ch_fasta
ch_dict
ch_fai

inputBamMutectTumorOnlyMode [bam,bai]

/*
================================================================================
=                      CREATE SOMATIC PANEL OF NORMALS                         =
================================================================================
*/

process run_mutect2_tumor_only_mode {

    tag "${bam.simpleName}"
    publishDir "${params.outdir}/MutectTumorOnlyMode", mode: 'copy'

    input:
        set file(bam),file(bai) from inputBamMutectTumorOnlyMode
        file(fasta) from ch_fasta
        file(dict) from ch_dict
        file(fastaFai) from ch_fai
        file(intervals) from ch_intervals

    output:
    file('*.vcf.gz') into inputVcfCreateGenomicsDB
    file('*.vcf.gz.tbi') into inputVcfIndexCreateGenomicsDB

    script:
    """
    gatk Mutect2 \
    -R ${ref} \
    -I ${bam} -normal ${bam.simpleName} \
    --max-mnp-distance 0 \
    -O ${bam.simpleName}.vcf.gz \
    -L $intervals \
    --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true'
    """
}


process create_GenomicsDB {

    tag ""
    publishDir "${params.outdir}/GenomicsDBImport", mode: 'copy'

    input:
        file("*") from inputVcfCreateGenomicsDB.collect()
        file("*") from inputVcfIndexCreateGenomicsDB.collect()
        file(fasta) from ch_fasta
        file(dict) from ch_dict
        file(fastaFai) from ch_fai
        file(intervals) from ch_intervals

    output:
    file("pon_db") into ponDbCreateSomaticPoN

    shell:
    //TODO Groovy > Bash
    '''
    echo -n "gatk GenomicsDBImport -R !{ref} --genomicsdb-workspace-path pon_db " > create_GenomicsDB.sh
    for vcf in $(ls *.vcf.gz); do
    echo -n "-V $vcf " >> create_GenomicsDB.sh
    done
    echo -n "-L !{intervals}" --merge-input-intervals --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' >> create_GenomicsDB.sh
    chmod ugo+xr create_GenomicsDB.sh
    bash create_GenomicsDB.sh
    '''
}

process create_somatic_PoN {
    
    tag "$germline_resource"
    publishDir "${params.outdir}/CreateSomaticPanelOfNormals", mode: 'copy'

    input:
        file(pon_db) from ponDbCreateSomaticPoN
        file(fasta) from ch_fasta
        file(dict) from ch_dict
        file(fastaFai) from ch_fai
        file(intervals) from ch_intervals
        file(germlineResource) from ch_germline_resource
        file(germlineResourceIndex) from ch_germline_resource_tbi

    output:
    set file("pon.vcf.gz"), file("pon.vcf.gz.tbi") into pon_vcf
    
    script:
    """
    gatk CreateSomaticPanelOfNormals \
    -R $ref \
    --germline-resource $germlineResource \
    -V gendb://$pon_db \
    -O pon.vcf.gz  \
    --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true'
    """
}

