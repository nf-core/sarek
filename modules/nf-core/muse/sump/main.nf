process MUSE_SUMP {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/83d1d3caa1b6ce54ce999e0061d7fe8acbe6788d5c7970574eff330ea819fb85/data'
        : 'community.wave.seqera.io/library/htslib_muse:9a4b9cb78c211f1e'}"

    input:
    tuple val(meta), path(muse_call_txt), path(ref_vcf), path(ref_vcf_tbi)

    output:
    tuple val(meta), path("*.vcf.gz"),     emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('muse'),  eval("MuSE --version | sed -e 's/MuSE, version //g' | sed -e 's/MuSE v//g'"), topic: versions, emit: versions_muse
    tuple val("${task.process}"), val('bgzip'), eval("bgzip --version | sed -n 's/bgzip (htslib) \\([0-9.]*\\)/\\1/p'"),      topic: versions, emit: versions_bgzip

    when:
    task.ext.when == null || task.ext.when

    script:
    // -G for WGS data and -E for WES data
    def args = task.ext.args ?: ''
    // args for bgzip
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // MuSE complains if the timestamp of the dbsnp VCF index is older than the timestamp of the VCF itself, so we need to touch it here
    """
    touch ${ref_vcf_tbi}

    MuSE \\
        sump \\
        ${args} \\
        -I ${muse_call_txt} \\
        -n ${task.cpus} \\
        -D ${ref_vcf} \\
        -O ${prefix}.vcf

    bgzip ${args2} --threads ${task.cpus} ${prefix}.vcf
    tabix -p vcf ${prefix}.vcf.gz
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    echo ${args2}
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
