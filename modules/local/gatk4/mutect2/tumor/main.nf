// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_MUTECT2_TUMOR {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.0.0--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.0.0--hdfd78af_1"
    }

    input:
    tuple val(meta), path(cram_normal), path(crai_normal), path(cram_tumor), path(crai_tumor), path(interval)
    path pon
    path ponIndex
    path dict
    path fasta
    path fai
    val no_intervals
    path(germline_resource)
    path(germline_resource_tbi)

    output:
    tuple val(meta), path("*.vcf"),       emit: vcf
    tuple val(meta), path("*.vcf.stats"), emit: vcf_stats
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def intervalsOptions = no_intervals ? "" : "-L ${interval}"
    def softClippedOption = params.ignore_soft_clipped_bases ? "--dont-use-soft-clipped-bases true" : ""
    def PON = params.pon ? "--panel-of-normals ${pon}" : ""
    if (!task.memory) {
        log.info '[GATK Mutect2] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    # Get raw calls
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        Mutect2 \
        -R ${fasta}\
        -I ${cramTumor} -tumor ${meta.tumor} \
        -I ${cramNormal} -normal ${meta.normal} \
        ${intervalsOptions} \
        ${softClippedOption} \
        --germline-resource ${germlineResource} \
        ${PON} \
        -O ${prefix}.vcf

    echo \$(gatk Mutect2 --version 2>&1) | sed 's/^.*(GATK) v//; s/ HTSJDK.*\$//' > ${software}.version.txt
    """
}
