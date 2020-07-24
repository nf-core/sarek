process SAMTOOLS_FAIDX {
    tag "${fasta}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        path fasta

    output:
        path "${fasta}.fai"

    script:
    """
    samtools faidx ${fasta}
    """
}
