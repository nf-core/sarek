process GATK4_GATHERPILEUPSUMMARIES {
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
    tuple val(meta), path(pileup)
    path  dict

    output:
    tuple val(meta), path("*.pileups.table"), emit: table
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_list = pileup.collect{ "--I $it" }.join(' ')

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK GatherPileupSummaries] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M" GatherPileupSummaries \\
        $input_list \\
        --O ${prefix}.pileups.table \\
        --sequence-dictionary $dict \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
