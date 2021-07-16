// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BWA_INDEX {
    tag "$fasta"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'index', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::bwa=0.7.17" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bwa:0.7.17--hed695b0_7"
    } else {
        container "quay.io/biocontainers/bwa:0.7.17--hed695b0_7"
    }

    input:
    path fasta

    output:
    path "bwa"          , emit: index
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    mkdir bwa
    bwa index $options.args $fasta -p bwa/${fasta.baseName}
    echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//' > ${software}.version.txt
    """
}
