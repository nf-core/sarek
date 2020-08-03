process SAMTOOLS_STATS {
    label 'cpus_2'

    tag "${meta.id}"

    publishDir "${params.outdir}/Reports/${meta.id}/SamToolsStats", mode: params.publish_dir_mode

    input:
        tuple val(meta), path(bam) 

    output:
        path ("${bam}.samtools.stats.out") 

    script:
    """
    samtools stats ${bam} > ${bam}.samtools.stats.out
    """
}