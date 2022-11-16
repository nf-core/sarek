process GATK4_BASERECALIBRATOR_SPARK {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.3.0.0 conda-forge::openjdk=8.0.312" : null)
    container 'broadinstitute/gatk:4.3.0.0'

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

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK BaseRecalibratorSpark] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" BaseRecalibratorSpark \\
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
