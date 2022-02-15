process GATK4_VARIANTRECALIBRATOR {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.5.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.5.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf) , path(tbi)
    path fasta
    path fai
    path dict
    tuple path(resvcfs), path(restbis), val(reslabels)

    output:
    tuple val(meta), path("*.recal")   , emit: recal
    tuple val(meta), path("*.idx")     , emit: idx
    tuple val(meta), path("*.tranches"), emit: tranches
    tuple val(meta), path("*plots.R")  , emit: plots, optional:true
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    refCommand = fasta ? "-R ${fasta} " : ''
    resourceCommand = '--resource:' + reslabels.join( ' --resource:')

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK VariantRecalibrator] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" VariantRecalibrator \\
        ${refCommand} \\
        -V ${vcf} \\
        -O ${prefix}.recal \\
        --tranches-file ${prefix}.tranches \\
        ${resourceCommand} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
