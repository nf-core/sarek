process BWAMEM2_INDEX {
    tag "${fasta}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/BWAIndex/${it}" : null }

    container "quay.io/biocontainers/bwa-mem2:2.0--he513fc3_1"

    conda (params.conda ? "bioconda::bwa-mem2=2.0" : null)

    input:
        path fasta
        val options

    output:
        path "${fasta}.*"

    script:
    """
    bwa-mem2 index ${fasta}
    """
}
