process GATK4_FILTERVARIANTTRANCHES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ce8f4142326abbb74e28b97633a200b09439c54480803343b1cf13dc19e51ae5/data'
        : 'community.wave.seqera.io/library/gatk4-lite:4.7.0.0--79918b8b2632f4b1'}"

    input:
    tuple val(meta), path(vcf), path(tbi), path(intervals)
    path resources
    path resources_index
    path fasta
    path fai
    path dict

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def resource_list = resources.collect { resource -> "--resource ${resource}" }.join(' ')

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK FilterVariantTranches] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        FilterVariantTranches \\
        --variant ${vcf} \\
        ${resource_list} \\
        --output ${prefix}.filtered.vcf.gz \\
        --tmp-dir . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip -c > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
