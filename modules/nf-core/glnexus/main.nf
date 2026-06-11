process GLNEXUS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6c/6cf504ad8e4ebda286609fa3c1a5f9af68dbca9ec06bb4428e219e84754bd140/data' :
        'community.wave.seqera.io/library/bcftools_glnexus:cf380f1a6410f606' }"

    input:
    tuple val(meta), path(gvcfs), path(custom_config)
    tuple val(meta2), path(bed)

    output:
    tuple val(meta), path("*.bcf"), emit: bcf
    tuple val("${task.process}"), val('glnexus'), eval("glnexus_cli 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+'"), topic: versions, emit: versions_glnexus

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions = bed ? "--bed ${bed}" : ""

    // Make list of GVCFs to merge
    def input = gvcfs.collect {vcf -> vcf.toString() }
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Glnexus] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    glnexus_cli \\
        --threads $task.cpus \\
        --mem-gbytes $avail_mem \\
        $regions \\
        $args \\
        ${input.join(' ')} \\
        > ${prefix}.bcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bcf
    """
}
