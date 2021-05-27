// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process INDEX_TARGET_BED {
    tag "$target_bed"
    label 'process_medium'
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:"false") }

    conda (params.enable_conda ? "bioconda::htslib=1.12" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/htslib:1.12--hd3b49d5_0"
    } else {
        container "quay.io/biocontainers/htslib:1.12--hd3b49d5_0"
    }

    input:
        path target_bed

    output:
        tuple path("${target_bed}.gz"), path("${target_bed}.gz.tbi")

    script:
    """
    bgzip --threads ${task.cpus} -c ${target_bed} > ${target_bed}.gz
    tabix ${target_bed}.gz
    """
}