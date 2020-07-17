process MERGE_BAM_MAPPED {
    label 'cpus'

    tag "${patient}-${sample}"

    input:
        tuple patient, sample, run, path(bam), path(bai)

    output:
        tuple patient, sample, path("${sample}.bam"), path("${sample}.bam.bai")

    script:
    """
    samtools merge --threads ${task.cpus} ${sample}.bam ${bam}
    samtools index ${sample}.bam
    """
}