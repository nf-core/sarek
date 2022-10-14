process SNAPALIGNER_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::snap-aligner=2.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snap-aligner:2.0.1--hd03093a_1':
        'quay.io/biocontainers/snap-aligner:2.0.1--hd03093a_1' }"

    input:
    tuple val(meta) , path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bai"), optional: true, emit: bai
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def subcmd = meta.single_end ? "single" : "paired"

    """
    INDEX=`dirname \$(find -L ./ -name "OverflowTable*")`
    [ -z "\$INDEX" ] && echo "Snap index files not found" 1>&2 && exit 1

    snap-aligner ${subcmd} \\
        \$INDEX \\
        ${reads.join(" ")} \\
        -o ${prefix}.bam \\
        -t ${task.cpus} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snapaligner: \$(snap-aligner 2>&1| head -n 1 | sed 's/^.*version //;s/.\$//')
    END_VERSIONS
    """
}
