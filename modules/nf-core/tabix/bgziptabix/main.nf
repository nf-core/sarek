process TABIX_BGZIPTABIX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92859404d861ae01afb87e2b789aebc71c0ab546397af890c7df74e4ee22c8dd/data' :
        'community.wave.seqera.io/library/htslib:1.21--ff8e28a189fbecaa' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.gz"), path("*.{tbi,csi}"), emit: gz_index
    tuple val("${task.process}"), val('tabix'), eval("tabix -h 2>&1 | grep -oP 'Version:\\s*\\K[^\\s]+'")   , topic: versions   , emit: versions_tabix
    tuple val("${task.process}"), val('bgzip'), eval("bgzip --version | sed '1!d;s/.* //'"), topic: versions, emit: versions_bgzip

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bgzip --threads ${task.cpus} -c $args $input > ${prefix}.${input.getExtension()}.gz
    tabix --threads ${task.cpus} $args2 ${prefix}.${input.getExtension()}.gz

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args2 = task.ext.args2 ?: ''
    def index = args2.contains("-C ") || args2.contains("--csi") ? "csi" : "tbi"
    """
    echo "" | gzip > ${prefix}.${input.getExtension()}.gz
    touch ${prefix}.${input.getExtension()}.gz.${index}

    """
}
