process MERGE_BAM {
    label 'cpus_8'

    tag "${meta.id}"

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("${meta.sample}.bam")

    //     when: !(params.no_intervals)

    script:
    """
    samtools merge --threads ${task.cpus} ${meta.sample}.bam ${bam}
    """
}