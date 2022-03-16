process GATK4_BASERECALIBRATOR_SPARK {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.3.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.3.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.3.0--hdfd78af_0' }"
    input:
    tuple val(meta), path(cram), path(crai), path(intervals_bed)
    path  fasta
    path  fasta_fai
    path  dict
    path  known_sites
    path  known_sites_tbi

    output:
    tuple val(meta), path("*.table"), emit: table
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args  ?: ''
    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK BaseRecalibratorSpark] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def intervals_command = intervals_bed ? "-L ${intervals_bed}" : ""
    def sites_command = known_sites.collect{"--known-sites ${it}"}.join(' ')
    """
    gatk BaseRecalibratorSpark  \
        -R $fasta \
        -I $cram \
        $sites_command \
        $intervals_command \
        --tmp-dir . \
        $args \
        -O ${prefix}.table \
        --spark-master local[${task.cpus}]

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
