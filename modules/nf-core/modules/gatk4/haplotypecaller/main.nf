process GATK4_HAPLOTYPECALLER {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.3.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.3.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.3.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(input), path(input_index)
    path fasta
    path fai
    path dict
    path dbsnp
    path dbsnp_tbi
    path interval

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.tbi")   , emit: tbi
    path "versions.yml"              , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def interval_option = interval ? "-L ${interval}" : ""
    def dbsnp_option    = dbsnp ? "-D ${dbsnp}" : ""
    def avail_mem       = 3
    if (!task.memory) {
        log.info '[GATK HaplotypeCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk \\
        --java-options "-Xmx${avail_mem}g" \\
        HaplotypeCaller \\
        -R $fasta \\
        -I $input \\
        ${dbsnp_option} \\
        ${interval_option} \\
        -O ${prefix}.vcf.gz \\
        $args \\
        --tmp-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
