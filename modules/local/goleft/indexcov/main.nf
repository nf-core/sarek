process GOLEFT_INDEXCOV {
label 'process_single'

container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/goleft:0.2.4--h9ee0642_1':
        'biocontainers/goleft:0.2.4--h9ee0642_1' }"

input:
        val(meta)
        path(bams)
        path(fasta)
        path(fai)
output:
	path("${prefix}/*"),emit:output
	path("versions.yml"),emit:versions
script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def input_files = bams.findAll{it.name.endsWith(".bam")} + bams.findAll{it.name.endsWith(".crai")}
    def extranormalize = input_files.any{it.name.endsWith(".crai")} ? " --extranormalize " : ""
"""
    goleft indexcov \\
        --fai "${fai}"  \\
        --directory "${prefix}" \\
        ${extranormalize} \\
        $args \\
        ${input_files.join(" ")}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        goleft: \$(goleft --version 2>&1 | head -n 1 | sed 's/^.*goleft Version: //')
    END_VERSIONS
"""

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir "${prefix}"
    touch "${prefix}/${prefix}-indexcov.bed.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        goleft: \$(goleft --version 2>&1 | head -n 1 | sed 's/^.*goleft Version: //')
    END_VERSIONS
    """
}




