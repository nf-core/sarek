process MERGE_BAM {
    label 'cpus_8'

    tag "${meta.id}"

    input:
        tuple val(meta), path(bam)//, path(bai) optional: true

    output:
        tuple val(meta), path("${meta.sample}.bam")//, path("${meta.sample}.bam.bai") optional: true

    //     when: !(params.no_intervals)
//    samtools merge --threads ${task.cpus} ${idSample}.bam ${bam}
//      samtools merge --threads ${task.cpus} ${idSample}.recal.bam ${bam}
    script:
    """
    samtools merge --threads ${task.cpus} ${meta.sample}.bam ${bam}
    """
}