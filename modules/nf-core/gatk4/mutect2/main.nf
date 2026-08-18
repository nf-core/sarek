process GATK4_MUTECT2 {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(input), path(input_index), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai), path(gzi)
    tuple val(meta4), path(dict)
    path alleles
    path alleles_tbi
    path germline_resource
    path germline_resource_tbi
    path panel_of_normals
    path panel_of_normals_tbi

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.stats"), emit: stats
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("${prefix}.f1r2.tar.gz"), emit: f1r2, optional: true
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def inputs = input.collect { vcf_ -> "--input ${vcf_}" }.join(" ")
    def interval_command = intervals ? "--intervals ${intervals}" : ""
    def pon_command = panel_of_normals ? "--panel-of-normals ${panel_of_normals}" : ""
    def gr_command = germline_resource ? "--germline-resource ${germline_resource}" : ""
    def a_command = alleles ? "--alleles ${alleles}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK Mutect2] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        Mutect2 \\
        ${inputs} \\
        --output ${prefix}.vcf.gz \\
        --reference ${fasta} \\
        ${pon_command} \\
        ${gr_command} \\
        ${a_command} \\
        ${interval_command} \\
        --tmp-dir . \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.vcf.gz.stats
    echo "" | gzip > ${prefix}.f1r2.tar.gz
    """
}
