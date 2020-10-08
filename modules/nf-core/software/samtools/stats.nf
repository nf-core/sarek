include { initOptions; saveFiles; getSoftwareName } from './../functions'

process SAMTOOLS_STATS {
    label 'cpus_2'

    tag "${meta.id}"

    publishDir "${params.outdir}/Reports/${meta.id}/SamToolsStats", mode: params.publish_dir_mode

    container "quay.io/biocontainers/samtools:1.10--h2e538c0_3"

    conda (params.conda ? "bioconda::samtools=1.10" : null)

    input:
        tuple val(meta), path(bam) 

    output:
        path ("${bam}.samtools.stats.out") 

    script:
    """
    samtools stats ${bam} > ${bam}.samtools.stats.out
    """
}