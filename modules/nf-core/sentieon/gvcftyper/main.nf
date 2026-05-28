process SENTIEON_GVCFTYPER {
    tag "${meta.id}"
    label 'process_high'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/73e9111552beb76e2ad3ad89eb75bed162d7c5b85b2433723ecb4fc96a02674a/data'
        : 'community.wave.seqera.io/library/sentieon:202503.02--def60555294d04fa'}"

    input:
    tuple val(meta), path(gvcfs), path(tbis), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dbsnp)
    tuple val(meta5), path(dbsnp_tbi)

    output:
    tuple val(meta), path("*.vcf.gz"),     emit: vcf_gz
    tuple val(meta), path("*.vcf.gz.tbi"), emit: vcf_gz_tbi
    tuple val("${task.process}"), val('sentieon'), eval('sentieon driver --version | sed "s/.*-//g"'), topic: versions, emit: versions_sentieon

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_genotyped"
    def gvcfs_input = '-v ' + gvcfs.join(' -v ')
    def dbsnp_cmd = dbsnp ? "--dbsnp ${dbsnp}" : ""
    def interval_command = intervals ? "--interval ${intervals}" : ""
    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64
        ? "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; "
        : ""
    """
    ${sentieonLicense}

    sentieon driver \\
        -r ${fasta} \\
        ${interval_command} \\
        --algo GVCFtyper \\
        ${gvcfs_input} \\
        ${dbsnp_cmd} \\
        ${prefix}.vcf.gz

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip >${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    """
}
