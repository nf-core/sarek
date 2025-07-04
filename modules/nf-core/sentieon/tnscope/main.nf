process SENTIEON_TNSCOPE {
    tag "$meta.id"
    label 'process_high'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/80/80ccb05eb4f1a193a3bd99c4da90f55f74ea6556c25f154e53e1ff5a6caa372d/data' :
        'community.wave.seqera.io/library/sentieon:202503--5e378058d837c58c' }"

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
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: index
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args     ?: ''
    def args2        = task.ext.args2    ?: ''
    def interval_str = intervals         ? "--interval ${intervals}" : ''
    def cosmic_str   = cosmic            ? "--cosmic ${cosmic}"      : ''
    def dbsnp_str    = dbsnp             ? "--dbsnp ${dbsnp}"        : ''
    def pon_str      = pon               ? "--pon ${pon}"            : ''
    def prefix       = task.ext.prefix   ?: "${meta.id}"
    def inputs       = input.collect{ "-i $it"}.join(" ")
    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64 ?
        "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; " :
        ""
    """
    $sentieonLicense


    sentieon driver \\
        -t $task.cpus \\
        -r $fasta \\
        $inputs \\
        $interval_str \\
        $args \\
        --algo TNscope \\
        $args2 \\
        $cosmic_str \\
        $dbsnp_str \\
        $pon_str \\
        ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g" )
    END_VERSIONS
    """
}
