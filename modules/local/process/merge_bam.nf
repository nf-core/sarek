include { initOptions; saveFiles; getSoftwareName } from './../../nf-core/software/functions'

params.options = [:]
def options    = initOptions(params.options)

environment = params.enable_conda ? "bioconda::samtools=1.10" : null
container = "quay.io/biocontainers/samtools:1.10--h2e538c0_3"
if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) container = "https://depot.galaxyproject.org/singularity/samtools:1.10--h2e538c0_3"

process MERGE_BAM {
    label 'cpus_8'

    tag "${meta.id}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda environment
    container container

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("${name}.bam"), emit: bam
        val meta,                            emit: tsv

    script:
    name = options.suffix ? "${meta.id}.${options.suffix}" : "${meta.id}"
    """
    samtools merge --threads ${task.cpus} ${name}.bam ${bam}
    """
}
