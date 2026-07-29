process CONTROLFREEC_MAKEGRAPH2 {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/control-freec:11.6b--hde5307d_3'
        : 'quay.io/biocontainers/control-freec:11.6b--hde5307d_3'}"

    input:
    tuple val(meta), path(ratio), path(baf)

    output:
    tuple val(meta), path("*_BAF.png"), emit: png_baf
    tuple val(meta), path("*_ratio.log2.png"), emit: png_ratio_log2
    tuple val(meta), path("*_ratio.png"), emit: png_ratio
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('controlfreec'), val("11.6b"), emit: versions_controlfreec, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def baf_cmd = baf ?: ""
    """
    cat \$(which makeGraph2.0.R) | R --slave --args ${args} ${ratio} ${baf_cmd}

    mv *_BAF.txt.png ${prefix}_BAF.png
    mv *_ratio.txt.log2.png ${prefix}_ratio.log2.png
    mv *_ratio.txt.png ${prefix}_ratio.png
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_BAF.png
    touch ${prefix}_ratio.log2.png
    touch ${prefix}_ratio.png
    """
}
