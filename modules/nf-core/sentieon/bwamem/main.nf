process SENTIEON_BWAMEM {
    tag "${meta.id}"
    label 'process_high'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/73e9111552beb76e2ad3ad89eb75bed162d7c5b85b2433723ecb4fc96a02674a/data'
        : 'community.wave.seqera.io/library/sentieon:202503.02--def60555294d04fa'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fasta_fai)

    output:
    tuple val(meta), path("${prefix}"), path("${prefix}.{bai,crai}"), emit: bam_and_bai
    tuple val("${task.process}"), val('bwa'), eval('sentieon bwa 2>&1 | sed -n "s/^Version: *//p"'), topic: versions, emit: versions_bwa
    tuple val("${task.process}"), val('sentieon'), eval('sentieon driver --version | sed "s/.*-//g"'), topic: versions, emit: versions_sentieon

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.bam"
    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64
        ? "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; "
        : ""

    """
    ${sentieonLicense}
    export bwt_max_mem="${(task.memory * 0.9).toGiga()}G"

    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

    sentieon bwa mem \\
        ${args} \\
        -t ${task.cpus} \\
        \$INDEX \\
        ${reads} \\
        | sentieon util sort -r ${fasta} -t ${task.cpus} -o ${prefix} --sam2bam -

    # Delete *.bai file if prefix ends with .cram
    if [[ "${prefix}" == *.cram ]]; then
        rm -f "${prefix}.bai"
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.bam"
    index = prefix.tokenize('.')[-1] == "bam" ? "bai" : "crai"
    """
    touch ${prefix}
    touch ${prefix}.${index}
    """
}
