process GATK4SPARK_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4-spark:4.4.0.0--hdfd78af_0':
        'biocontainers/gatk4-spark:4.4.0.0--hdfd78af_0' }"

    if (workflow.containerEngine == 'docker') {
        containerOptions = '-u $(id -u):$(id -g) -v /etc/passwd:/etc/passwd:ro'
    }

    input:
    tuple val(meta), path(bam)
    path  fasta
    path  fasta_fai
    path  dict

    output:
    tuple val(meta), path("${prefix}"),     emit: output
    tuple val(meta), path("${prefix}.bai"), emit: bam_index, optional:true
    tuple val(meta), path("*.metrics"),     emit: metrics, optional: true
    path "versions.yml"               ,     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def input_list = bam.collect{"--input $it"}.join(' ')

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK MarkDuplicatesSpark] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        MarkDuplicatesSpark \\
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
