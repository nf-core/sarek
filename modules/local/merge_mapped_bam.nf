process MERGE_BAM_MAPPED {
    label 'cpus'

    tag "${idPatient}-${idSample}"

    input:
        tuple idPatient, idSample, idRun, path(bam), path(bai)// from multiple

    output:
        tuple idPatient, idSample, path("${idSample}.bam") //into bam_mapped_merged

    script:
    """
    samtools merge --threads ${task.cpus} ${idSample}.bam ${bam}
    """
}