process SENTIEON_TNSCOPE {
    tag "${meta.id}"
    label 'process_high'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/73e9111552beb76e2ad3ad89eb75bed162d7c5b85b2433723ecb4fc96a02674a/data'
        : 'community.wave.seqera.io/library/sentieon:202503.02--def60555294d04fa'}"

    input:
    tuple val(meta), path(input), path(input_index), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dbsnp)
    tuple val(meta5), path(dbsnp_tbi)
    tuple val(meta6), path(pon)
    tuple val(meta7), path(pon_tbi)
    tuple val(meta8), path(cosmic)
    tuple val(meta9), path(cosmic_tbi)

    output:
    tuple val(meta), path("*.vcf.gz"),     emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: index
    tuple val("${task.process}"), val('sentieon'), eval('sentieon driver --version | sed "s/.*-//g"'), topic: versions, emit: versions_sentieon

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def interval_str = intervals ? "--interval ${intervals}" : ''
    def cosmic_str = cosmic ? "--cosmic ${cosmic}" : ''
    def dbsnp_str = dbsnp ? "--dbsnp ${dbsnp}" : ''
    def pon_str = pon ? "--pon ${pon}" : ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def inputs = input.collect {in -> "-i ${in}" }.join(" ")
    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64
        ? "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; "
        : ""
    """
    ${sentieonLicense}


    sentieon driver \\
        -t ${task.cpus} \\
        -r ${fasta} \\
        ${inputs} \\
        ${interval_str} \\
        ${args} \\
        --algo TNscope \\
        ${args2} \\
        ${cosmic_str} \\
        ${dbsnp_str} \\
        ${pon_str} \\
        ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
