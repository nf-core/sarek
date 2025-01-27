process GATK4SPARK_BASERECALIBRATOR {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4-spark:4.5.0.0--hdfd78af_0':
        'biocontainers/gatk4-spark:4.5.0.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(input), path(input_index), path(intervals)
    path  fasta
    path  fai
    path  dict
    path  known_sites
    path  known_sites_tbi

    output:
    tuple val(meta), path("*.table"), emit: table
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def interval_command = intervals ? "--intervals $intervals" : ""
    def sites_command = known_sites.collect{"--known-sites $it"}.join(' ')

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK BaseRecalibratorSpark] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        BaseRecalibratorSpark \\
        --input $input \\
        --output ${prefix}.table \\
        --reference $fasta \\
        $interval_command \\
        $sites_command \\
        --spark-master local[${task.cpus}] \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
