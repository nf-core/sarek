process SAMTOOLS_MERGE_CRAM {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.14" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0' :
        'quay.io/biocontainers/samtools:1.14--hb421002_0' }"

    input:
    tuple val(meta), path(crams)
    path  fasta

    output:
    tuple val(meta), path("*.cram"), emit: cram
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    """
    samtools merge -@${task.cpus} --reference ${fasta} ${prefix}.cram $crams

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
