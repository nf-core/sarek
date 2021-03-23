include { initOptions; saveFiles; getSoftwareName } from './../functions'

params.options = [:]
def options    = initOptions(params.options)

process MSISENSOR_SCAN {
    tag "$fasta"
    label 'process_high'
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:"false") }

    conda (params.enable_conda ? "bioconda::msisensor=0.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/msisensor:0.5--hb3646a4_2"
    } else {
        container "quay.io/biocontainers/msisensor:0.5--hb3646a4_2"
    }

    input:
        path fasta
        path fai

    output:
        path "microsatellites.list"

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "msisensor_${ioptions.suffix}" : "msisensor_"
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "$ioptions.args" variable
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    """
    msisensor scan -d ${fasta} -o microsatellites.list
    """
}