process BWA_INDEX {
    tag "$fasta"
    // NOTE requires 5.37N memory where N is the size of the database
    // source: https://bio-bwa.sourceforge.net/bwa.shtml#8
    memory { 6.B * fasta.size() }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bf/bf7890f8d4e38a7586581cb7fa13401b7af1582f21d94eef969df4cea852b6da/data' :
        'community.wave.seqera.io/library/bwa_htslib_samtools:56c9f8d5201889a4' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bwa")  , emit: index
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    def args   = task.ext.args ?: ''
    """
    mkdir bwa
    bwa \\
        index \\
        $args \\
        -p bwa/${prefix} \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    mkdir bwa

    touch bwa/${prefix}.amb
    touch bwa/${prefix}.ann
    touch bwa/${prefix}.bwt
    touch bwa/${prefix}.pac
    touch bwa/${prefix}.sa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}
