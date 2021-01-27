include { initOptions; saveFiles; getSoftwareName } from './../functions'

params.options = [:]
def options    = initOptions(params.options)

environment = params.enable_conda ? "bioconda::msisensor=0.5" : null
container = "quay.io/biocontainers/msisensor:0.5--hb3646a4_2"
if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) container = "https://depot.galaxyproject.org/singularity/msisensor:0.5--hb3646a4_2"

process MSISENSOR_MSI {
    tag "${meta.id}"
    
    label 'CPUS_1'
    label 'MEMORY_MAX'
    
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda environment
    container container

    input:
        tuple val(meta), path(bam_normal), path(bai_normal), path(bam_tumor), path(bai_tumor)
        path msisensor_scan

    output:
        tuple val(meta), path("*.list")

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "msisensor_${meta.id}${ioptions.suffix}" : "msisensor_${meta.id}"
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "$ioptions.args" variable
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    """
    msisensor msi -d ${msisensor_scan} \
                  -b 4 \
                  -t ${bam_tumor} \
                  -n ${bam_normal} \
                  -o ${prefix}

    mv ${prefix}          ${prefix}.list
    mv ${prefix}_dis      ${prefix}_dis.list
    mv ${prefix}_germline ${prefix}_germline.list
    mv ${prefix}_somatic  ${prefix}_somatic.list
    """
}