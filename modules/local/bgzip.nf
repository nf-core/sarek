process BGZIP {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::pbgzip=2016.08.04" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbgzip:2016.08.04--h36cd882_2':
        'quay.io/biocontainers/pbgzip:2016.08.04--h36cd882_2' }"

    input:
    tuple val(meta), path(vcf_gz)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pbgzip -d $vcf_gz -n${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bgzip: \$(echo \$(pbgzip --version 2>&1) | sed 's/^.*pbgzip //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
