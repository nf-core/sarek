process RTGTOOLS_FORMAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rtg-tools:3.12.1--hdfd78af_0':
        'biocontainers/rtg-tools:3.12.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(input1), path(input2), path(sam_rg)

    output:
    tuple val(meta), path("*.sdf"), emit: sdf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def single = meta.containsKey("single_end") ? meta.single_end : true

    def input = single ? "${input1}" : "--left ${input1} --right ${input2}"
    def rg = sam_rg ? "--sam-rg ${sam_rg}" : ""

    def avail_mem = "3G"
    if (!task.memory) {
        log.info '[RTG format] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue() + "M"
    }

    """
    rtg RTG_MEM=${avail_mem} format \\
        ${args} \\
        ${rg} \\
        --output ${prefix}.sdf \\
        ${input}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rtg-tools: \$(echo \$(rtg version | head -n 1 | awk '{print \$4}'))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = "3G"
    if (!task.memory) {
        log.info '[RTG format] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue() + "M"
    }
    """
    touch ${prefix}.sdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rtg-tools: \$(echo \$(rtg version | head -n 1 | awk '{print \$4}'))
    END_VERSIONS
    """
}
