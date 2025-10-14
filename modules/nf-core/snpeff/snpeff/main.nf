process SNPEFF_SNPEFF {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/37a5be8afdc2e112fd1bb4aa226e009a437e6537f50a51ed909fe2018fd62e99/data' :
        'community.wave.seqera.io/library/snpeff:5.3.0a--ca8e0b08f227a463' }"

    input:
    tuple val(meta),  path(vcf)
    val db
    tuple val(meta2), path(cache)

    output:
    tuple val(meta), path("*.ann.vcf"),     emit: vcf
    tuple val(meta), path("*.csv"),         emit: report
    tuple val(meta), path("*.html"),        emit: summary_html
    tuple val(meta), path("*.genes.txt"),   emit: genes_txt
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def avail_mem = 6144
    if (!task.memory) {
        log.info('[snpEff] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cache_command = cache ? "-dataDir \${PWD}/${cache}" : ""
    """
    snpEff \\
        -Xmx${avail_mem}M \\
        ${db} \\
        ${args} \\
        -csvStats ${prefix}.csv \\
        ${cache_command} \\
        ${vcf} \\
        > ${prefix}.ann.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ann.vcf
    touch ${prefix}.csv
    touch ${prefix}.html
    touch ${prefix}.genes.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """
}
