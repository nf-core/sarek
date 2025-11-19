process RBT_VCFSPLIT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rust-bio-tools:0.42.2--h4458251_1':
        'biocontainers/rust-bio-tools:0.42.2--h4458251_1' }"

    input:
    tuple val(meta), path(vcf)
    val(numchunks)

    output:
    tuple val(meta), path("*.bcf"), emit: bcfchunks
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunks = numchunks ? (numchunks - 1) : 15
    """
    rbt vcf-split \\
        ${vcf} \\
        ${prefix}.{0..${chunks}}.bcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Rust-Bio-Tools: \$(rbt --version | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunks = numchunks ? (numchunks - 1) : 15
    def bcf_files = (0..chunks).collect { "${prefix}.${it}.bcf" }.join(' ')
    """
    touch ${bcf_files}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Rust-Bio-Tools: \$(rbt --version | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/')
    END_VERSIONS
    """
}
