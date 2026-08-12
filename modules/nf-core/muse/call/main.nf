process MUSE_CALL {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9f/9f0ebb574ef5eed2a6e034f1b2feea6c252d1ab0c8bc5135a669059aa1f4d2ca/data'
        : 'community.wave.seqera.io/library/muse:6637291dcbb0bdb8'}"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(reference)

    output:
    tuple val(meta), path("*.MuSE.txt"), emit: txt
    tuple val("${task.process}"), val('muse'),  eval("MuSE --version | sed -e 's/MuSE, version //g' | sed -e 's/MuSE v//g'"), topic: versions, emit: versions_muse

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    MuSE \\
        call \\
        ${args} \\
        -f ${reference} \\
        -O ${prefix}  \\
        -n ${task.cpus} \\
        ${tumor_bam}    \\
        ${normal_bam}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    touch ${prefix}.MuSE.txt
    """
}
