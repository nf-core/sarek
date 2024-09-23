process DEEPVARIANT_MAKEEXAMPLES {
    tag "$meta.id"
    label 'process_high'

    //Conda is not supported at the moment
    container "nf-core/deepvariant:1.6.1"

    input:
    tuple val(meta), path(input), path(index), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(gzi)
    tuple val(meta5), path(par_bed)

    output:
    tuple val(meta), path("${prefix}.examples.tfrecord-*-of-*.gz{,.example_info.json}"),    emit: examples
    tuple val(meta), path("${prefix}.gvcf.tfrecord-*-of-*.gz"),        emit: gvcf
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DEEPVARIANT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def regions = intervals ? "--regions ${intervals}" : ""
    def par_regions = par_bed ? "--par_regions_bed=${par_bed}" : ""

    """
    seq 0 ${task.cpus - 1} | parallel -q --halt 2 --line-buffer /opt/deepvariant/bin/make_examples \\
        --mode calling \\
        --ref "${fasta}" \\
        --reads "${input}" \\
        --examples "./${prefix}.examples.tfrecord@${task.cpus}.gz" \\
        --gvcf "./${prefix}.gvcf.tfrecord@${task.cpus}.gz" \\
        ${regions} \\
        ${par_regions} \\
        ${args} \\
        --task {}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant_makeexamples: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf -v SHARD_COUNT "%04d" ${task.cpus}
    for i in \$( seq -f "%04g" 0 ${task.cpus-1} )
    do
        touch ${prefix}.examples.tfrecord-\$i-of-\$SHARD_COUNT.tfrecord.gz{,.example_info.json}
        touch ${prefix}.gvcf.tfrecord-\$i-of-\$SHARD_COUNT.tfrecord.gz
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant_makeexamples: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    """
}
