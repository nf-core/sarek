process GATK4_CALCULATECONTAMINATION {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(pileup), path(matched)

    output:
    tuple val(meta), path('*.contamination.table'), emit: contamination
    tuple val(meta), path('*.segmentation.table'), emit: segmentation, optional: true
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def matched_command = matched ? "--matched-normal ${matched}" : ''

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK CalculateContamination] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        CalculateContamination \\
        --input ${pileup} \\
        --output ${prefix}.contamination.table \\
        ${matched_command} \\
        --tmp-dir . \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.contamination.table
    touch ${prefix}.segmentation.table
    """
}
