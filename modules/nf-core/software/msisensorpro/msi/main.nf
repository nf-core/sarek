// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MSISENSORPRO_MSI {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::msisensor-pro=1.1.a" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/msisensor-pro:1.1.a--hb3646a4_0"
    } else {
        container "quay.io/biocontainers/msisensor-pro:1.1.a--hb3646a4_0"
    }

    input:
        tuple val(meta), path(bam_normal), path(bai_normal), path(bam_tumor), path(bai_tumor)
        path msisensorpro_scan

    output:
        tuple val(meta), path("msisensorpro_*.list")

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}.${options.suffix}" : "${meta.id}"
    """
    msisensor-pro msi \\
        -d $msisensorpro_scan \\
        -n $bam_normal \\
        -t $bam_tumor \\
        -o $prefix \\
        -b $task.cpus \\
        $options.args

    mv ${prefix}          msisensorpro_${prefix}.list
    mv ${prefix}_dis      msisensorpro_${prefix}_dis.list
    mv ${prefix}_germline msisensorpro_${prefix}_germline.list
    mv ${prefix}_somatic  msisensorpro_${prefix}_somatic.list
    """
}