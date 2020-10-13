include { initOptions; saveFiles; getSoftwareName } from './../../nf-core/software/functions'

environment = params.conda ? "bioconda::samtools=1.10" : null
container = "quay.io/biocontainers/samtools:1.10--h2e538c0_3"
if (workflow.containerEngine == 'singularity') container = "https://depot.galaxyproject.org/singularity/samtools:1.10--h2e538c0_3"

process MERGE_BAM {
    label 'cpus_8'

    tag "${meta.id}"

    conda environment
    container container

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
