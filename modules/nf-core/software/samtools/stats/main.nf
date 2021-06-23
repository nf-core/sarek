// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SAMTOOLS_STATS {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::samtools=1.12" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.12--hd5e65b6_0"
    } else {
        container "quay.io/biocontainers/samtools:1.12--hd5e65b6_0"
    }

    input:
    tuple val(meta), path(cram), path(crai)
    path(reference)

    output:
    tuple val(meta), path("*.stats"), emit: stats
    path  "*.version.txt"           , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    export REF_PATH="$(pwd)/ref/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s"
    export REF_CACHE="$(pwd)/ref/cache/%2s/%2s/%s"
    samtools stats --ref-seq ${reference} $cram > ${cram}.stats
    echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' > ${software}.version.txt
    """
}
