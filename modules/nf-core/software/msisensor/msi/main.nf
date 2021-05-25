// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MSISENSOR_MSI {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::msisensor=0.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/msisensor:0.5--hb3646a4_2"
    } else {
        container "quay.io/biocontainers/msisensor:0.5--hb3646a4_2"
    }

    input:
        tuple val(meta), path(bam_normal), path(bai_normal), path(bam_tumor), path(bai_tumor)
        path msisensor_scan

    output:
        tuple val(meta), path("*.list")

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}.${options.suffix}" : "${meta.id}"
    """
    msisensor msi \
        -d ${msisensor_scan} \
        -b 4 \
        -t ${bam_tumor} \
        -n ${bam_normal} \
        -o ${prefix} \
        $options.args \

    mv ${prefix}          msisensor_${prefix}.list
    mv ${prefix}_dis      msisensor_${prefix}_dis.list
    mv ${prefix}_germline msisensor_${prefix}_germline.list
    mv ${prefix}_somatic  msisensor_${prefix}_somatic.list
    """
}