process MUSE_SUMP {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f9/f951cafbcec988c70d08bec486d40d902db4114bb1aaee89fd68155fdc3273da/data':
        'community.wave.seqera.io/library/muse_tabix:cf4c19969fa62f63' }"

    input:
    tuple val(meta) , path(muse_call_txt)
    tuple val(meta2), path(ref_vcf), path(ref_vcf_tbi)

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: '' // -G for WGS data and -E for WES data
    def args2  = task.ext.args2  ?: '' // args for bgzip
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    MuSE \\
        sump \\
        $args \\
        -I $muse_call_txt \\
        -n $task.cpus    \\
        -D $ref_vcf \\
        -O ${prefix}.vcf

    bgzip $args2 --threads $task.cpus ${prefix}.vcf
    tabix -p vcf ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MuSE: \$( MuSE --version | sed -e "s/MuSE, version //g" | sed -e "s/MuSE v//g")
        bgzip: \$( bgzip --version | sed -n 's/bgzip (htslib) \\([0-9.]*\\)/\\1/p' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MuSE: \$( MuSE --version | sed -e "s/MuSE, version //g" | sed -e "s/MuSE v//g")
        bgzip: \$( bgzip --version | sed -n 's/bgzip (htslib) \\([0-9.]*\\)/\\1/p' )
    END_VERSIONS
    """
}
