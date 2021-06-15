// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_MARKDUPLICATES_SPARK {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gatk4-spark=4.2.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4-spark:4.2.0.0--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/gatk4-spark:4.2.0.0--hdfd78af_1"
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
    def metrics  = use_metrics ? "-M ${prefix}.bam.metrics" :''
    def bams     = bam.collect(){ x -> "-I ".concat(x.toString()) }.join(" ")

    """
    gatk MarkDuplicatesSpark \\
        ${bams} \\
        $metrics \
        --tmp-dir . \\
        --create-output-bam-index true \\
        --spark-master local[${task.cpus}] \\
        -O ${prefix}.bam \\
        $options.args

    echo \$(gatk MarkDuplicatesSpark --version 2>&1) | sed 's/^.*(GATK) v//; s/ HTSJDK.*\$//' > ${software}.version.txt
    """
}