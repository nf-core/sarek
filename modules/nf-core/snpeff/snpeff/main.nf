process SNPEFF_SNPEFF {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/30/30669e5208952f30d59d0d559928772f082830d01a140a853fff13a2283a17b0/data'
        : 'community.wave.seqera.io/library/snpeff:5.4.0a--eaf6ce30125b2b17'}"

    input:
    tuple val(meta), path(vcf)
    val db
    tuple val(meta2), path(cache)

    output:
    tuple val(meta), path("*.ann.vcf"), emit: vcf
    tuple val(meta), val("${task.process}"), val('snpeff'), path("*.csv"), topic: multiqc_files, emit: report
    tuple val(meta), val("${task.process}"), val('snpeff'), path("*.html"), topic: multiqc_files, emit: summary_html
    tuple val(meta), val("${task.process}"), val('snpeff'), path("*.genes.txt"), topic: multiqc_files, emit: genes_txt
    tuple val("${task.process}"), val('snpeff'), eval("snpEff -version 2>&1 | cut -f 2 -d '\t'"), topic: versions, emit: versions_snpeff

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
        -XX:-UsePerfData \\
        ${db} \\
        ${args} \\
        -csvStats ${prefix}.csv \\
        ${cache_command} \\
        ${vcf} \\
        > ${prefix}.ann.vcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ann.vcf
    touch ${prefix}.csv
    touch ${prefix}.html
    touch ${prefix}.genes.txt
    """
}
