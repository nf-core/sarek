process GATK4_MUTECT2_MERGE {
    //tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.3.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.3.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.3.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(cram_normal), path(crai_normal), path(cram_tumor), path(crai_tumor), path(interval)
    path  fasta_fai
    path  fasta
    path  dict
    path  pon
    path  pon_tbi
    path  germline_resource
    path  germline_resource_tbi

    output:
    tuple val(meta), path("*.vcf")      , emit: vcf
    tuple val(meta), path("*.vcf.stats"), emit: vcf_stats
    path "versions.yml"                 , emit: versions

    script:
    def args = task.ext.args  ?: ''
    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK Mutect2] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def intervals_command = intervals_bed ? "-L ${intervals_bed}" : ""
    def pon_command = pon ? "--panel-of-normals ${pon}" : ""
    // def softClippedOption = params.ignore_soft_clipped_bases ? "--dont-use-soft-clipped-bases true" : ""
    if (!task.memory) {
        log.info '[GATK Mutect2] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    // def softClippedOption = params.ignore_soft_clipped_bases ? "--dont-use-soft-clipped-bases true" : ""
    """
    # Get raw calls
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        Mutect2 \
        -R $fasta \
        -I $cram_tumor -tumor ${meta.tumor} \
        -I $cram_normal -normal ${meta.normal} \
        $intervals_command \
        $pon_command \
        $args \
        --germline-resource $germline_resource \
        -O ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
