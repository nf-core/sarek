process SENTIEON_APPLYVARCAL {
    tag "$meta.id"
    label 'process_low'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/sentieon:202308.02--ffce1b7074ce9924' :
        'nf-core/sentieon:202308.02--c641bc397cbf79d5' }"

    input:
    tuple val(meta), path(vcf), path(vcf_tbi), path(recal), path(recal_index), path(tranches)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.tbi")   , emit: tbi
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64 ?
        "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; " :
        ""
    """
    $sentieonLicense

    sentieon driver -r ${fasta}  --algo ApplyVarCal \\
        -v $vcf \\
        --recal $recal \\
        --tranches_file $tranches \\
        $args \\
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
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
