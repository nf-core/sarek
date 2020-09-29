process SAMTOOLS_FAIDX {
    tag "${fasta}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    container "quay.io/biocontainers/samtools:1.10--h2e538c0_3"

    conda (params.conda ? "bioconda::samtools=1.10" : null)

    input:
        path fasta

    output:
        path "${fasta}.fai"

    script:
    """
    samtools faidx ${fasta}
    """
}
