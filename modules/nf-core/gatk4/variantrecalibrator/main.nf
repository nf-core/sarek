process GATK4_VARIANTRECALIBRATOR {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::gatk4=4.4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0':
        'biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi) // input vcf and tbi of variants to recalibrate
    path resource_vcf   // resource vcf
    path resource_tbi   // resource tbi
    val labels          // string (or list of strings) containing dedicated resource labels already formatted with '--resource:' tag
    path  fasta
    path  fai
    path  dict

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
    def reference_command = fasta ? "--reference $fasta " : ''
    def labels_command = labels.join(' ')

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK VariantRecalibrator] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        VariantRecalibrator \\
        --variant $vcf \\
        --output ${prefix}.recal \\
        --tranches-file ${prefix}.tranches \\
        $reference_command \\
        --tmp-dir . \\
        $labels_command \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.recal
    touch ${prefix}.idx
    touch ${prefix}.tranches
    touch ${prefix}plots.R

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
