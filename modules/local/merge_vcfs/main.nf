process MERGE_VCFS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.5.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.5.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)
    path  fasta_fai
    path  target_bed

    output:
    tuple val(meta), path("${prefix}.vcf.gz")    , emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: tbi
    path  "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args  ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def target_options   = target_bed ? "-t ${target_bed}" : ""

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK ApplyBQSRSpark] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    ls *.vcf.gz > ${prefix}.vcf.list
    gatk --java-options "-Xmx${avail_mem}g" MergeVcfs -I ${prefix}.vcf.list -O ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
