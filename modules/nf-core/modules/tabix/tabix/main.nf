// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process TABIX_TABIX {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::tabix=1.11' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/tabix:1.11--hdfd78af_0"
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
