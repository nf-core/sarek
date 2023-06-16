process GATK4_APPLYBQSR {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::gatk4=4.4.0.0"
    container { NfcoreTemplate.getContainer(
        docker: 'biocontainers/gatk4:4.4.0.0--py36hdfd78af_0',
        singularity: 'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0',
        registry: 'quay.io',
        use_full_uri: task.ext.container_full_uri ?: false
    )}

    input:
    tuple val(meta), path(input), path(input_index), path(bqsr_table), path(intervals)
    path  fasta
    path  fai
    path  dict

    output:
    tuple val(meta), path("*.bam") , emit: bam,  optional: true
    tuple val(meta), path("*.cram"), emit: cram, optional: true
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def interval_command = intervals ? "--intervals $intervals" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK ApplyBQSR] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M" ApplyBQSR \\
        --input $input \\
        --output ${prefix}.${input.getExtension()} \\
        --reference $fasta \\
        --bqsr-recal-file $bqsr_table \\
        $interval_command \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
