process GATK4_CALCULATECONTAMINATION {
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
    tuple val(meta), path(pileup), path(matched)

    output:
    tuple val(meta), path('*.contamination.table'), emit: contamination
    tuple val(meta), path('*.segmentation.table') , emit: segmentation, optional:true
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def matched_command = matched ? "--matched-normal $matched" : ''

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK CalculateContamination] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M" CalculateContamination \\
        --input $pileup \\
        --output ${prefix}.contamination.table \\
        $matched_command \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
