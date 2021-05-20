// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_HAPLOTYPECALLER {
    tag "$meta.id"
    label 'process_medium'
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
    tuple val(meta), path(bam), path(bai), path(interval)
    path dbsnp
    tuple val(meta_dbsnp), path(dbsnp_tbi)
    path dict
    path fasta
    path fai
    val no_intervals

    output:
    tuple val(meta), path("*.vcf.gz"),                 emit: vcf
    tuple val(meta), path(interval), path("*.vcf.gz"), emit: interval_vcf
    tuple val(meta), path("*.tbi"),                    emit: tbi
    path "*.version.txt",                              emit: version

    script:
    def software  = getSoftwareName(task.process)
    def prefix    = options.suffix ? "${interval.baseName}_${meta.id}${options.suffix}" : "${interval.baseName}_${meta.id}"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK HaplotypeCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    def intervalsOptions = no_intervals ? "" : "-L ${interval}"
    def dbsnpOptions = params.dbsnp ? "--D ${dbsnp}" : ""
    """
    gatk \\
        --java-options "-Xmx${avail_mem}g" \\
        HaplotypeCaller \\
        -R $fasta \\
        -I $bam \\
        ${intervalsOptions} \\
        ${dbsnpOptions} \\
        -O ${prefix}.vcf.gz \\
        $options.args

    gatk --version | grep Picard | sed "s/Picard Version: //g" > ${software}.version.txt
    """
}
