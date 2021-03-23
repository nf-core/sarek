// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process TABIX_TABIX {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::tabix=0.2.6" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/tabix:0.2.6--ha92aebf_0"
    } else {
        container "quay.io/biocontainers/tabix:0.2.6--ha92aebf_0"
    }

    input:
    tuple val(meta), path(tab)

    output:
    tuple val(meta), path("*.tbi"), emit: tbi
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    tabix $options.args $tab

    echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/(.*\$//' > ${software}.version.txt
    """
}
