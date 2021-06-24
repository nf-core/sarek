// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.0.0--0"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.0.0--0"
    }

    input:
    tuple val(meta), path(bam)//, path(bai)
    val use_metrics

    output:
    tuple val(meta), path("*.bam"), path("*.bai"),       emit: bam
    tuple val(meta), path("*.metrics"), optional : true, emit: metrics
    path "*.version.txt",                                emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def metrics  = use_metrics ? "M=${prefix}.metrics" :''
    def bams     = bam.collect(){ x -> "INPUT=".concat(x.toString()) }.join(" ")
    if (!task.memory) {
        log.info '[GATK MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk MarkDuplicates \\
        ${bams} \\
        $metrics \
        TMP_DIR=. \\
        ASSUME_SORT_ORDER=coordinate \\
        CREATE_INDEX=true \\
        O=${prefix}.bam \\
        $options.args

    mv ${prefix}.bai ${prefix}.bam.bai

    echo \$(gatk MarkDuplicates --version 2>&1) | sed 's/^.*(GATK) v//; s/ HTSJDK.*\$//' > ${software}.version.txt
    """
}