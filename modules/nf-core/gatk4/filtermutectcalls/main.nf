process GATK4_FILTERMUTECTCALLS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.5.0.0--py36hdfd78af_0':
        'biocontainers/gatk4:4.5.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(vcf_tbi), path(stats), path(orientationbias), path(segmentation), path(table), val(estimate)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)

    output:
    tuple val(meta), path("*.vcf.gz")            , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi")        , emit: tbi
    tuple val(meta), path("*.filteringStats.tsv"), emit: stats
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def orientationbias_command = orientationbias ? orientationbias.collect{"--orientation-bias-artifact-priors $it"}.join(' ') : ''
    def segmentation_command    = segmentation    ? segmentation.collect{"--tumor-segmentation $it"}.join(' ')                  : ''
    def estimate_command        = estimate        ? " --contamination-estimate ${estimate} "                                    : ''
    def table_command           = table           ? table.collect{"--contamination-table $it"}.join(' ')                        : ''

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK FilterMutectCalls] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        FilterMutectCalls \\
        --variant $vcf \\
        --output ${prefix}.vcf.gz \\
        --reference $fasta \\
        $orientationbias_command \\
        $segmentation_command \\
        $estimate_command \\
        $table_command \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.vcf.gz.filteringStats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
