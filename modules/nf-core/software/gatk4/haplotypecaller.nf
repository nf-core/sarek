include { initOptions; saveFiles; getSoftwareName } from './../functions'

params.options = [:]
def options    = initOptions(params.options)

process GATK4_HAPLOTYPECALLER {
    tag "$meta.id"
    label 'process_medium'
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.0.0--0"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.0.0--0"
    }

    input:
    tuple val(meta), path(bam), path(bai), file(interval)
    path dbsnp
    tuple val(meta_dbsnp), path(dbsnp_tbi)
    path dict
    path fasta
    path fai
    val no_intervals

    output:
    tuple val(meta), path("${interval.baseName}_${meta.id}.g.vcf"),                 emit: gvcf
    tuple val(meta), path(interval), path("${interval.baseName}_${meta.id}.g.vcf"), emit: interval_gvcf
    path "*.version.txt" , emit: version

    script:
    def software = getSoftwareName(task.process)
    intervalsOptions = no_intervals ? "" : "-L ${interval}"
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

    echo \$(gatk HaplotypeCaller --version 2>&1) | sed 's/^.*(GATK) v//; s/ HTSJDK.*\$//' > ${software}.version.txt
    """
}
