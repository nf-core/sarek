process GATK4_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::gatk4=4.4.0.0 bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d9e7bad0f7fbc8f4458d5c3ab7ffaaf0235b59fb:f857e2d6cc88d35580d01cf39e0959a68b83c1d9-0':
        'biocontainers/mulled-v2-d9e7bad0f7fbc8f4458d5c3ab7ffaaf0235b59fb:f857e2d6cc88d35580d01cf39e0959a68b83c1d9-0' }"

    input:
    tuple val(meta), path(bam)
    path  fasta
    path  fasta_fai

    output:
    tuple val(meta), path("*cram"),     emit: cram,  optional: true
    tuple val(meta), path("*bam"),      emit: bam,   optional: true
    tuple val(meta), path("*.crai"),    emit: crai,  optional: true
    tuple val(meta), path("*.bai"),     emit: bai,   optional: true
    tuple val(meta), path("*.metrics"), emit: metrics
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.bam"

    // If the extension is CRAM, then change it to BAM
    prefix_bam = prefix.tokenize('.')[-1] == 'cram' ? "${prefix.substring(0, prefix.lastIndexOf('.'))}.bam" : prefix

    def input_list = bam.collect{"--INPUT $it"}.join(' ')
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    // Using samtools and not Markduplicates to compress to CRAM speeds up computation:
    // https://medium.com/@acarroll.dna/looking-at-trade-offs-in-compression-levels-for-genomics-tools-eec2834e8b94
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        MarkDuplicates \\
        $input_list \\
        --OUTPUT ${prefix_bam} \\
        --METRICS_FILE ${prefix}.metrics \\
        --TMP_DIR . \\
        ${reference} \\
        $args

    # If cram files are wished as output, the run samtools for conversion
    if [[ ${prefix} == *.cram ]]; then
        samtools view -Ch -T ${fasta} -o ${prefix} ${prefix_bam}
        rm ${prefix_bam}
        samtools index ${prefix}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
