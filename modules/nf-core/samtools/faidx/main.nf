process SAMTOOLS_FAIDX {
    tag "$fasta"
    label 'process_single'

    conda "bioconda::samtools=1.17"
    container { NfcoreTemplate.getContainer(
        docker: 'biocontainers/samtools:1.17--h00cdaf9_0',
        singularity: 'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0',
        registry: 'quay.io',
        use_full_uri: task.ext.container_full_uri ?: false
    )}

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(fai)

    output:
    tuple val(meta), path ("*.{fa,fasta}") , emit: fa , optional: true
    tuple val(meta), path ("*.fai")        , emit: fai, optional: true
    tuple val(meta), path ("*.gzi")        , emit: gzi, optional: true
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools \\
        faidx \\
        $fasta \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def match = (task.ext.args =~ /-o(?:utput)?\s(.*)\s?/).findAll()
    def fastacmd = match[0] ? "touch ${match[0][1]}" : ''
    """
    ${fastacmd}
    touch ${fasta}.fai

    cat <<-END_VERSIONS > versions.yml

    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
