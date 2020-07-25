process MERGE_BAM {
    label 'cpus_8'

    tag "${meta.id}"

    input:
        tuple val(meta), path(bam), path(bai)

    output:
        tuple val(meta), path("${meta.sample}.bam"), path("${meta.sample}.bam.bai")

    script:
    """
    samtools merge --threads ${task.cpus} ${meta.sample}.bam ${bam}
    samtools index ${meta.sample}.bam
    """
}