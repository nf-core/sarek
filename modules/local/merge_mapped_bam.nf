process MergeBamMapped {
    label 'cpus'

    tag "${idPatient}-${idSample}"

    input:
        tuple idPatient, idSample, idRun, path(bam) // from multiple

    output:
        tuple idPatient, idSample, path("${idSample}.bam") //into bam_mapped_mer

    script:
    """
    samtools merge --threads ${task.cpus} ${idSample}.bam ${bam}
    """
}