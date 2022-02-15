process GATK4_APPLYBQSR {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.5.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.5.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(input), path(input_index), path(bqsr_table), path(intervals)
    path  fasta
    path  fai
    path  dict

    output:
    tuple val(meta), path("*.bam"),  emit: bam, optional: true
    tuple val(meta), path("*.cram"), emit: cram, optional: true
    path "versions.yml"           ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def interval = intervals ? "-L ${intervals}" : ""
    def file_type = input.getExtension()

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK ApplyBQSR] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" ApplyBQSR \\
        -R $fasta \\
        -I $input \\
        --bqsr-recal-file $bqsr_table \\
        $interval \\
        --tmp-dir . \\
        -O ${prefix}.${file_type} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
