process GATK4SPARK_BASERECALIBRATOR {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/49/498aea9c9bcaf736b9fb2a01366c1b7b38ccc0d38143178afc325d6a93241447/data'
        : 'community.wave.seqera.io/library/gatk4-spark:4.6.2.0--8b5cd67ee60a714e'}"

    input:
    tuple val(meta), path(input), path(input_index), path(intervals)
    path fasta
    path fai
    path dict
    path known_sites
    path known_sites_tbi

    output:
    tuple val(meta), path("*.table"), emit: table
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def interval_command = intervals ? "--intervals ${intervals}" : ""
    def sites_command = known_sites.collect { vcf -> "--known-sites ${vcf}" }.join(' ')

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK BaseRecalibratorSpark] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        BaseRecalibratorSpark \\
        --input ${input} \\
        --output ${prefix}.table \\
        --reference ${fasta} \\
        ${interval_command} \\
        ${sites_command} \\
        --spark-master local[${task.cpus}] \\
        --tmp-dir . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.table
    """
}
