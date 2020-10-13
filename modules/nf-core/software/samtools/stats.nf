include { initOptions; saveFiles; getSoftwareName } from './../functions'

environment = params.conda ? "bioconda::samtools=1.10" : null
container = "quay.io/biocontainers/samtools:1.10--h2e538c0_3"
if (workflow.containerEngine == 'singularity') container = "https://depot.galaxyproject.org/singularity/samtools:1.10--h2e538c0_3"

process SAMTOOLS_STATS {
    label 'cpus_2'

    tag "${meta.id}"

    publishDir "${params.outdir}/Reports/${meta.id}/SamToolsStats", mode: params.publish_dir_mode

    conda environment
    container container

    input:
        tuple val(meta), path(bam) 

    output:
        path ("${bam}.samtools.stats.out") 

    script:
    """
    samtools stats ${bam} > ${bam}.samtools.stats.out
    """
}