include { initOptions; saveFiles; getSoftwareName } from './../functions'

params.options = [:]
def options    = initOptions(params.options)

process GATK_HAPLOTYPECALLER {
    label 'MEMORY_SINGLECPU_TASK_SQ'
    label 'CPUS_2'

    tag "${meta.id}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::gatk4=4.1.9.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.1.9.0--py39_0"
    } else {
        container "quay.io/biocontainers/gatk4:4.1.9.0--py39_0"
    }

    input:
        tuple val(meta), path(bam), path(bai), file(interval)
        path dbsnp
        path dbsnpIndex
        path dict
        path fasta
        path fai

    output:
        tuple val(meta), path("${interval.baseName}_${meta.id}.g.vcf"),                 emit: gvcf
        tuple val(meta), path(interval), path("${interval.baseName}_${meta.id}.g.vcf"), emit: interval_gvcf


    script:
    intervalsOptions = params.no_intervals ? "" : "-L ${interval}"
    dbsnpOptions = params.dbsnp ? "--D ${dbsnp}" : ""
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        ${intervalsOptions} \
        ${dbsnpOptions} \
        -O ${interval.baseName}_${meta.id}.g.vcf \
        -ERC GVCF
    """
}
