process BCFTOOLS_STATS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(regions)
    tuple val(meta3), path(targets)
    tuple val(meta4), path(samples)
    tuple val(meta5), path(exons)
    tuple val(meta6), path(fasta)

    output:
    tuple val(meta), path("*stats.txt"), emit: stats
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions_file = regions ? "--regions-file ${regions}" : ""
    def targets_file = targets ? "--targets-file ${targets}" : ""
    def samples_file = samples ? "--samples-file ${samples}" : ""
    def reference_fasta = fasta ? "--fasta-ref ${fasta}" : ""
    def exons_file = exons ? "--exons ${exons}" : ""
    """
    bcftools stats \\
        ${args} \\
        ${regions_file} \\
        ${targets_file} \\
        ${samples_file} \\
        ${reference_fasta} \\
        ${exons_file} \\
        ${vcf} > ${prefix}.bcftools_stats.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bcftools_stats.txt
    """
}
