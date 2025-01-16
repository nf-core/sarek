process EXPANSIONHUNTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/expansionhunter:5.0.0--hf366f20_0' :
        'biocontainers/expansionhunter:5.0.0--hf366f20_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    tuple val(meta4), path(variant_catalog)

    output:
    tuple val(meta), path("*.vcf.gz")        , emit: vcf
    tuple val(meta), path("*.json.gz")       , emit: json
    tuple val(meta), path("*_realigned.bam") , emit: bam
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    ExpansionHunter \\
        ${args} \\
        --reads ${bam} \\
        --output-prefix ${prefix} \\
        --reference ${fasta} \\
        --variant-catalog ${variant_catalog}

    bgzip --threads ${task.cpus} ${args2} ${prefix}.vcf
    bgzip --threads ${task.cpus} ${args2} ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        expansionhunter: \$( echo \$(ExpansionHunter --version 2>&1) | head -1 | sed 's/^.*ExpansionHunter v//')
        bgzip: \$(echo \$(bgzip -h 2>&1) | sed 's/^.*Version: //;s/Usage:.*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    echo "" | gzip > ${prefix}.json.gz
    touch ${prefix}_realigned.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        expansionhunter: \$( echo \$(ExpansionHunter --version 2>&1) | head -1 | sed 's/^.*ExpansionHunter v//')
        bgzip: \$(echo \$(bgzip -h 2>&1) | sed 's/^.*Version: //;s/Usage:.*//')
    END_VERSIONS
    """
}
