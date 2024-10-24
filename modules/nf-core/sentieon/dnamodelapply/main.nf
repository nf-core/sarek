process SENTIEON_DNAMODELAPPLY {
    tag "$meta.id"
    label 'process_high'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a6/a64461f38d76bebea8e21441079e76e663e1168b0c59dafee6ee58440ad8c8ac/data' :
        'community.wave.seqera.io/library/sentieon:202308.03--59589f002351c221' }"

    input:
    tuple val(meta), path(vcf), path(idx)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(ml_model)

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: index
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64 ?
        "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; " :
        ""
    """
    $sentieonLicense

    sentieon driver \\
        -t $task.cpus \\
        -r $fasta \\
        $args \\
        --algo DNAModelApply \\
        --model $ml_model \\
        -v $vcf \\
        ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g" )
    END_VERSIONS
    """
}
