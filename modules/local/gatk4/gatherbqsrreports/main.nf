process GATK4_GATHERBQSRREPORTS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.3.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.3.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.3.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(recal)

    output:
    tuple val(meta), path("*.recal.table"), emit: table
    path "*.recal.table"                  , emit: report
    path "versions.yml"                   , emit: versions

    script:
    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK GatherBQSRReports] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def input = recal.collect{"-I ${it}"}.join(' ')
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        GatherBQSRReports \
        ${input} \
        --tmp-dir . \
        -O ${meta.sample}.recal.table

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
