process GATK4_MARKDUPLICATES_SPARK {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::gatk4=4.3.0.0 conda-forge::openjdk=8.0.312" : null)
    container 'broadinstitute/gatk:4.3.0.0'

    input:
    tuple val(meta), path(bam)
    path  fasta
    path  fasta_fai
    path  dict

    output:
    tuple val(meta), path("${prefix}"),     emit: output
    tuple val(meta), path("*.metrics"),     emit: metrics, optional: true
    path "versions.yml"               ,     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def input_list = bam.collect{"--input $it"}.join(' ')


    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK MarkDuplicatesSpark] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" MarkDuplicatesSpark \\
        $input_list \\
        --output $prefix \\
        --reference $fasta \\
        --spark-master local[${task.cpus}] \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        openjdk: \$(echo \$(java -version 2>&1) | grep version | sed 's/\"//g' | cut -f3 -d ' ')
    END_VERSIONS
    """
}
