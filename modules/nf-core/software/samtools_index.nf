process SAMTOOLS_INDEX {
   label 'cpus_8'

    tag "${meta.id}"

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path(bam), path("*.bai")
//    samtools index ${idSample}.bam

    //    samtools index ${idSample}.recal.bam

    script:
    """
    samtools index $bam
    """
}