process SAMTOOLS_SORT_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.14" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "versions.yml"          , emit: versions

    // tuple val(meta), path("*sorted.bam"), path("*.bai")  , optional:true, emit: bam_bai
    // tuple val(meta), path("*sorted.bam"), path("*.csi")  , optional:true, emit: bam_csi
    // path  "versions.yml"          , emit: versions

    script:
    """
    samtools sort -@ $task.cpus -o ${meta.id}.sorted.bam -T $meta.id $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
