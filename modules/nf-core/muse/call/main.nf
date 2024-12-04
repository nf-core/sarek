process MUSE_CALL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9f/9f0ebb574ef5eed2a6e034f1b2feea6c252d1ab0c8bc5135a669059aa1f4d2ca/data':
        'community.wave.seqera.io/library/muse:6637291dcbb0bdb8' }"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("*.MuSE.txt"), emit: txt
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    MuSE \\
        call \\
        $args \\
        -f $reference \\
        -O ${prefix}  \\
        -n $task.cpus \\
        $tumor_bam    \\
        $normal_bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MuSE: \$( MuSE --version | sed -e "s/MuSE, version //g" )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.MuSE.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MuSE: \$( MuSE --version | sed -e "s/MuSE, version //g" )
    END_VERSIONS
    """
}
