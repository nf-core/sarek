process XENGSORT_INDEX {
    tag "$host_fasta"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/xengsort:2.2.1--pyhdfd78af_0'
        : 'quay.io/biocontainers/xengsort:2.2.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(host_fasta, stageAs: "host/*")
    tuple val(meta2), path(graft_fasta, stageAs: "graft/*")
    val nobjects

    output:
    tuple val(meta2), path("*.hash")          , emit: hash
    tuple val(meta2), path("*.info")          , emit: info
    tuple val("${task.process}"), val('xengsort'), eval("xengsort --version"), topic: versions, emit: versions_xengsort

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cpus = task.cpus as int

    // The presented thread allocation behavior is taken from the tools source code
    // https://gitlab.com/genomeinformatics/xengsort/-/blob/master/xengsort/xengsort/xengsort_index.py?ref_type=heads
    def subtables = Math.max([(cpus / 2) as int - 1, cpus - 3, 19].min(), 1)
    if ((subtables % 2) == 0) {
        subtables += 1
    }

    def read_threads = Math.ceil(subtables / 10) as int

    def split_threads = 2 * read_threads

    if ((subtables + read_threads + split_threads) >= cpus) {
        read_threads = 1
        split_threads = 2
    }
    """
    xengsort index \\
        $args \\
        --threads-split ${split_threads} \\
        --threads-read ${read_threads} \\
        --subtables ${subtables} \\
        --index ${prefix} \\
        --host ${host_fasta} \\
        --graft ${graft_fasta} \\
        --nobjects ${nobjects}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.info
    touch ${prefix}.hash
    """
}
