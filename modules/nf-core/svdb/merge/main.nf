process SVDB_MERGE {
    tag "$meta.id"
    label 'process_medium'
    conda (params.enable_conda ? "bioconda::svdb=2.6.1 bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c8daa8f9d69d3c5a1a4ff08283a166c18edb0000:56d0a468970fbb474d92f0591abcf677757fb370-0':
        'quay.io/biocontainers/mulled-v2-c8daa8f9d69d3c5a1a4ff08283a166c18edb0000:56d0a468970fbb474d92f0591abcf677757fb370-0' }"

    input:
    tuple val(meta), path(vcfs)
    val (priority)

    output:
    tuple val(meta), path("*_sv_merge.vcf.gz"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input  = "${vcfs.join(" ")}"
    def prio   = ""
    if(priority) {
        prio = "--priority ${priority.join(',')}"
        input = ""
        for (int index = 0; index < vcfs.size(); index++) {
            input += " ${vcfs[index]}:${priority[index]}"
        }
    }
    """
    svdb \\
        --merge \\
        $args \\
        $prio \\
        --vcf $input \\
        > ${prefix}_sv_merge.vcf
    bgzip ${prefix}_sv_merge.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_sv_merge.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
