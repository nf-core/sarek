process IRONQC {
    tag "$meta.id"
    label 'process_medium'

    // Container built by Wave from co-located Dockerfile (no container directive needed)

    input:
    tuple val(meta), path(cram), path(crai), path(bed)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.stats")                       , emit: stats
    tuple val(meta), path("*.mosdepth.global.dist.txt")    , emit: global_dist
    tuple val(meta), path("*.mosdepth.region.dist.txt")    , optional: true, emit: region_dist
    tuple val(meta), path("*.mosdepth.summary.txt")        , emit: summary
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args   ?: ""
    """
    ironqc bundle \\
        ${cram} \\
        --reference ${fasta} \\
        --fai ${fai} \\
        --threads ${task.cpus} \\
        --prefix ${prefix} \\
        --indexcov-dir indexcov \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ironqc: \$(ironqc --version 2>&1 | head -1 | sed 's/ironqc //')
    END_VERSIONS
    """
}
