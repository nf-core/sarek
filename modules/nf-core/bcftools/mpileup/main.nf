process BCFTOOLS_MPILEUP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    tuple val(meta), path(bam), path(intervals_mpileup, stageAs: 'mpileup_intervals/*'), path(intervals_call, stageAs: 'call_intervals/*')
    tuple val(meta2), path(fasta), path(fai)
    val save_mpileup

    output:
    tuple val(meta), path("*vcf.gz"), emit: vcf
    tuple val(meta), path("*.{tbi,csi}"), emit: index, optional: true
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
    def intervals_mpileup_cmd = intervals_mpileup ? "-T ${intervals_mpileup}" : ""
    def intervals_call_cmd = intervals_call ? "-T ${intervals_call}" : ""
    """
    echo "${meta.id}" > sample_name.list

    bcftools \\
        mpileup \\
        --fasta-ref ${fasta} \\
        ${args} \\
        ${bam} \\
        ${intervals_mpileup_cmd} \\
        ${mpileup} \\
        | bcftools call --output-type v ${args2} ${intervals_call_cmd} \\
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
