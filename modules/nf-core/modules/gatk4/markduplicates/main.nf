process GATK4_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.5.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.5.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*.bam")    , emit: bam
    tuple val(meta), path("*.bai")    , emit: bai
    tuple val(meta), path("*.metrics"), emit: metrics
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bam_list = bams.collect(){ bam -> "--INPUT ".concat(bam.toString()) }.join(" ")
    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" MarkDuplicates \\
        $bam_list \\
        --METRICS_FILE ${prefix}.metrics \\
        --TMP_DIR . \\
        --CREATE_INDEX true \\
        --OUTPUT ${prefix}.bam \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
