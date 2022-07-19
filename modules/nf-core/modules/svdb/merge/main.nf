process SVDB_MERGE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::svdb=2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svdb:2.6.1--py39h5371cbf_0':
        'quay.io/biocontainers/svdb:2.6.1--py39h5371cbf_0' }"

    input:
    tuple val(meta), path(vcfs)
    val (priority)

    output:
    tuple val(meta), path("*_sv_merge.vcf"), emit: vcf
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

    cp ${prefix}_sv_merge.vcf /home/centos/git/sarek/outputsave/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_sv_merge.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
    END_VERSIONS
    """
}
