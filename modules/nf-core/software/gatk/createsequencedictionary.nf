process GATK_CREATESEQUENCEDICTIONARY {
    tag "${fasta}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    container "quay.io/biocontainers/gatk4-spark:4.1.8.1--0"

    conda (params.conda ? "bioconda::gatk4-spark=4.1.8.1" : null)

    input:
        path fasta

    output:
        path "${fasta.baseName}.dict"

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        CreateSequenceDictionary \
        --REFERENCE ${fasta} \
        --OUTPUT ${fasta.baseName}.dict
    """
}