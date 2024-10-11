process SAMTOOLS_COLLATEFASTQ {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(fasta)
    val(interleave)

    output:
    tuple val(meta), path("*_{1,2}.fq.gz")     , optional:true, emit: fastq
    tuple val(meta), path("*_interleaved.fq")  , optional:true, emit: fastq_interleaved
    tuple val(meta), path("*_other.fq.gz")     , emit: fastq_other
    tuple val(meta), path("*_singleton.fq.gz") , optional:true, emit: fastq_singleton
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""
    def output =    (interleave && ! meta.single_end) ? "> ${prefix}_interleaved.fq"                     :
                    meta.single_end                   ? "-1 ${prefix}_1.fq.gz -s ${prefix}_singleton.fq.gz" :
                    "-1 ${prefix}_1.fq.gz -2 ${prefix}_2.fq.gz -s ${prefix}_singleton.fq.gz"

    """
    samtools collate \\
        $args \\
        --threads $task.cpus \\
        ${reference} \\
        -O \\
        $input \\
        . |

    samtools fastq \\
        $args2 \\
        --threads $task.cpus \\
        ${reference} \\
        -0 ${prefix}_other.fq.gz \\
        $output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def empty = "echo '' | gzip "
    def singletoncommand = "${empty}> ${prefix}_singleton.fq.gz"
    def interleavecommand = interleave && !meta.single_end ? "${empty}> ${prefix}_interleaved.fq.gz" : ""
    def output1command = !interleave ? "${empty}> ${prefix}_1.fq.gz" : ""
    def output2command = !interleave && !meta.single_end ? "${empty}> ${prefix}_2.fq.gz" : ""

    """
    ${output1command}
    ${output2command}
    ${interleavecommand}
    ${singletoncommand}
    ${empty}> ${prefix}_other.fq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
