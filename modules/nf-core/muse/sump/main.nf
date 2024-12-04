process MUSE_SUMP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/35/3567f6162ff718c648175c5e7b5f848eaa27811d0cb3ad53def8f0a1c8893efa/data':
        'community.wave.seqera.io/library/muse_tabix:df58ca78bd9447b7' }"

    input:
    tuple val(meta), path(muse_call_txt)
    tuple val(meta2), path(ref_vcf), path(ref_vcf_tbi)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MuSE: \$( MuSE --version | sed -e "s/MuSE, version //g" )
        bgzip: \$( bgzip --version | sed -e "s/bgzip (htslib) //g" )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MuSE: \$( MuSE --version | sed -e "s/MuSE, version //g" )
        bgzip: \$( bgzip --version | sed -e "s/bgzip (htslib) //g" )
    END_VERSIONS
    """
}
