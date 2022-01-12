process GATK4_MARKDUPLICATES_SPARK {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.4.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.4.1--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.4.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    path  fasta
    path  fasta_fai
    path  dict
    val   format //either "bam" or "cram"

    output:
    tuple val(meta), path("*.${format}"), emit: output
    path "versions.yml"                 , emit: versions

    script:
    def args = task.ext.args  ?: ''
    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK MarkDuplicatesSpark] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def bams = bam.collect(){ x -> "-I ".concat(x.toString()) }.join(" ")
    """
    gatk  \
        MarkDuplicatesSpark \
        $bams \
        -O ${prefix}.${format} \
        --reference ${fasta} \
        --tmp-dir . \
        --spark-master local[${task.cpus}] \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
