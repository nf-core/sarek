process MERGE_BAM {
    label 'cpus_8'

    tag "${meta.id}"
    //TODO publishDir
    
    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("${meta.sample}.bam")

    script:
    """
    samtools merge --threads ${task.cpus} ${meta.sample}.bam ${bam}
    """
    //TODO Naming?
    //samtools merge --threads ${task.cpus} ${idSample}.recal.bam ${bam}
}