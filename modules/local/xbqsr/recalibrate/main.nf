process XBQSR_RECALIBRATE {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        "docker://ghcr.io/ikabadzhov/xbqsr:${task.ext.gpu_vendor ?: (params.gpu_vendor ?: 'latest')}" :
        "ghcr.io/ikabadzhov/xbqsr:${task.ext.gpu_vendor ?: (params.gpu_vendor ?: 'latest')}" }"

    input:
    tuple val(meta), path(input), path(input_index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(known_sites)

    output:
    tuple val(meta), path("*.table"), emit: table
    tuple val("${task.process}"), val('bqsr'), eval("bqsr --version 2>&1 | sed 's/bqsr //' || echo 1.0.0"), topic: versions, emit: versions_bqsr

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("X-BQSR module does not support Conda. Please use Docker / Singularity / Apptainer with -profile gpu instead.")
    }
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sites  = known_sites instanceof List ? known_sites.join(',') : "${known_sites}"
    """
    bqsr recalibrate \\
        --bam-in ${input} \\
        --recal-out ${prefix}.table \\
        --ref ${fasta} \\
        --known-sites ${sites} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.table
    """
}
