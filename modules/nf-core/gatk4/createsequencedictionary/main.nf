process GATK4_CREATESEQUENCEDICTIONARY {
    tag "$fasta"
    label 'process_medium'

    conda "bioconda::gatk4=4.4.0.0"
    container { NfcoreTemplate.getContainer(
        docker: 'biocontainers/gatk4:4.4.0.0--py36hdfd78af_0',
        singularity: 'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0',
        registry: 'quay.io',
        use_full_uri: task.ext.container_full_uri ?: false
    )}

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*.dict')  , emit: dict
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def avail_mem = 6144
    if (!task.memory) {
        log.info '[GATK CreateSequenceDictionary] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M" CreateSequenceDictionary \\
        --REFERENCE $fasta \\
        --URI $fasta \\
        --TMP_DIR . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta.baseName}.dict

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
