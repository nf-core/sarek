// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './../../nf-core/software/functions'

params.options = [:]
def options    = initOptions(params.options)

process SPLIT_FASTQ{
    label 'process_medium'
    //publishDir "${params.outdir}",
    //    mode: params.publish_dir_mode,
    //    saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'split_reads_seqkit', publish_id:'') }
    
    conda (params.enable_conda ? "bioconda::seqkit=0.13.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/seqkit:0.13.2--0"
    } else {
        container "quay.io/biocontainers/seqkit:0.13.2--0"
    }

    input:
        tuple val(name), path(read1), path(read2)

    output:
        tuple val(name), path ("*.gz")

    script:
    def software = getSoftwareName(task.process)
    """
    seqkit split2 --threads ${task.cpus} -1 $read1 -2 $read2 $options.args
    echo \$(seqkit --version 2>&1) > ${software}.version.txt
    """
}