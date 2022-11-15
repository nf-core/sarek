process SAMTOOLS_COLLATEFASTQ {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.16.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(fasta)
    val(interleave)

    output:
    tuple val(meta), path("*_{1,2}.fq.gz")          , optional:true, emit: fastq
    tuple val(meta), path("*_interleaved.fq.gz")    , optional:true, emit: fastq_interleaved
    tuple val(meta), path("*_other.fq.gz")          , emit: fastq_other
    tuple val(meta), path("*_singleton.fq.gz")      , optional:true, emit: fastq_singleton
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""
    def output =    (interleave && ! meta.single_end) ? "> ${prefix}_interleaved.fq.gz"                     :
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
}
