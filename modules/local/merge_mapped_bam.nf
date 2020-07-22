process MERGE_BAM_MAPPED {
    label 'cpus'

    tag "${patient}-${sample}"

    input:
        tuple val(patient), val(sample), val(run), path(bam), path(bai)

    output:
        tuple val(patient), val(sample), path("${sample}.bam"), path("${sample}.bam.bai")

    script:
    """
    samtools merge --threads ${task.cpus} ${sample}.bam ${bam}
    samtools index ${sample}.bam
    """
}