process GATK4_BASERECALIBRATOR {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.5.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.5.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(input), path(input_index), path(intervals)
    path fasta
    path fai
    path dict
    path knownSites
    path knownSites_tbi

    output:
    tuple val(meta), path("*.table"), emit: table
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def intervalsCommand = intervals ? "-L ${intervals}" : ""
    def sitesCommand = knownSites.collect{"--known-sites ${it}"}.join(' ')

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK BaseRecalibrator] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    """
    gatk --java-options "-Xmx${avail_mem}g" BaseRecalibrator  \
        -R $fasta \
        -I $input \
        $sitesCommand \
        $intervalsCommand \
        --tmp-dir . \
        $args \
        -O ${prefix}.table

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
