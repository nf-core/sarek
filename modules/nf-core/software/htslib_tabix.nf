// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process HTSLIB_TABIX {
    tag "${vcf}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:"false") }

    conda (params.enable_conda ? "bioconda::tabix=0.2.6" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/tabix:0.2.6--ha92aebf_0"
    } else {
        container "quay.io/biocontainers/tabix:0.2.6--ha92aebf_0"
    }

    input:
        path vcf

    output:
        path "${vcf}.tbi"

    script:
    def software = getSoftwareName(task.process)
    """
    tabix -p vcf ${vcf}

    echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/(.*\$//' > ${software}.version.txt
    """
}
