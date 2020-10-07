process MERGE_BAM {
    label 'cpus_8'

    tag "${meta.id}"

    container "quay.io/biocontainers/samtools:1.10--h2e538c0_3"

    conda (params.conda ? "bioconda::samtools=1.10" : null)

    input:
        tuple val(meta), path(bam)
        val options

    output:
        tuple val(meta), path("${name}.bam"), emit: bam
        val meta,                            emit: tsv

    script:
    name = options.suffix ? "${meta.id}.${options.suffix}" : "${meta.id}"
    """
    samtools merge --threads ${task.cpus} ${name}.bam ${bam}
    """
}
