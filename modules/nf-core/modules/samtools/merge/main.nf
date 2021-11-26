process SAMTOOLS_MERGE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.14" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0' :
        'quay.io/biocontainers/samtools:1.14--hb421002_0' }"

    input:
    tuple val(meta), path(input_files)
    path fasta

    output:
    tuple val(meta), path("${prefix}.bam"),  optional:true, emit: bam
    tuple val(meta), path("${prefix}.cram"), optional:true, emit: cram
    path  "versions.yml"                                  , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def file_type = input_files[0].getExtension()
    def reference = fasta ? "--reference ${fasta}" : ""
    """
    samtools merge --threads ${task.cpus-1} $args ${reference} ${prefix}.${file_type} $input_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
