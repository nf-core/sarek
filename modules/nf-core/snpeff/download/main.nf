process SNPEFF_DOWNLOAD {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/30/30669e5208952f30d59d0d559928772f082830d01a140a853fff13a2283a17b0/data'
        : 'community.wave.seqera.io/library/snpeff:5.4.0a--eaf6ce30125b2b17'}"

    input:
    tuple val(meta), val(snpeff_db)

    output:
    tuple val(meta), path('snpeff_cache'), emit: cache
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
    """
    snpEff \\
        -Xmx${avail_mem}M \\
        download ${snpeff_db} \\
        -dataDir \${PWD}/snpeff_cache \\
        ${args}
    """

    stub:
    """
    mkdir -p snpeff_cache/${snpeff_db}

    touch snpeff_cache/${snpeff_db}/sequence.I.bin
    touch snpeff_cache/${snpeff_db}/sequence.bin
    """
}
