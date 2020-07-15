process BWAMEM2_INDEX {
    tag "${fasta}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/BWAIndex/${it}" : null }

    input:
        path file(fasta)

    output:
        path file("${fasta}.*")

    script:
    """
    bwa-mem2 index ${fasta}
    """
}
