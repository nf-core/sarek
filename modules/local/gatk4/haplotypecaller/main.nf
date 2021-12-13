process GATK4_HAPLOTYPECALLER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.3.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.3.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.3.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(cram), path(crai), path(intervals_bed)
    path  fasta
    path  fasta_fai
    path  dict
    path  dbsnp
    path  dbsnp_tbi

    output:
    tuple val(meta), path("*.vcf")                , emit: vcf
    tuple val(meta), path("*.vcf"), path(intervals_bed), emit: interval_vcf
    path "versions.yml"                           , emit: versions

    script:
    def args = task.ext.args  ?: ''
    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK HaplotypeCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def intervals_command = intervals_bed ? "-L ${intervals_bed}" : ""
    def sites_command = dbsnp ? "--D ${dbsnp}" : ""
    //TODO allow ploidy argument here since we allow it for the cnv callers? or is this covered with options? Might unintuitive to use
    """
    gatk \\
        --java-options "-Xmx${avail_mem}g" \\
        HaplotypeCaller \\
        -R $fasta \\
        -I $cram \\
        ${sites_command} \\
        ${intervals_command} \\
        -O ${prefix}.vcf \\
        --tmp-dir . \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
