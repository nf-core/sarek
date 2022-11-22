process GATK4_LEARNREADORIENTATIONMODEL {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.3.0.0--py36hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(f1r2)

    output:
    tuple val(meta), path("*.tar.gz"), emit: artifactprior
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_list = f1r2.collect{"--input $it"}.join(' ')

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK LearnReadOrientationModel] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" LearnReadOrientationModel \\
        $input_list \\
        --output ${prefix}.tar.gz \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
