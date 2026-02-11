process BCFTOOLS_MPILEUP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f'}"

    input:
    tuple val(meta), path(bam), path(intervals)
    tuple val(meta2), path(fasta)
    val save_mpileup

    output:
    tuple val(meta), path("*vcf.gz"), emit: vcf
    tuple val(meta), path("*vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("*stats.txt"), emit: stats
    tuple val(meta), path("*.mpileup.gz"), emit: mpileup, optional: true
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mpileup = save_mpileup ? "| tee ${prefix}.mpileup" : ""
    def bgzip_mpileup = save_mpileup ? "bgzip ${prefix}.mpileup" : ""
    def intervals_cmd = intervals ? "-T ${intervals}" : ""
    """
    echo "${meta.id}" > sample_name.list

    bcftools \\
        mpileup \\
        --fasta-ref ${fasta} \\
        ${args} \\
        ${bam} \\
        ${intervals_cmd} \\
        ${mpileup} \\
        | bcftools call --output-type v ${args2} \\
        | bcftools reheader --samples sample_name.list \\
        | bcftools view --output-file ${prefix}.vcf.gz --output-type z ${args3}

    ${bgzip_mpileup}

    tabix -p vcf -f ${prefix}.vcf.gz

    bcftools stats ${prefix}.vcf.gz > ${prefix}.bcftools_stats.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bcftools_stats.txt
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    echo "" | gzip > ${prefix}.mpileup.gz
    """
}
