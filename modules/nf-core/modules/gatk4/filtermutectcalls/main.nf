process GATK4_FILTERMUTECTCALLS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.5.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.5.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(stats), path(orientationbias), path(segmentation), path(contaminationfile), val(contaminationest)
    path fasta
    path fai
    path dict

    output:
    tuple val(meta), path("*.vcf.gz")            , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi")        , emit: tbi
    tuple val(meta), path("*.filteringStats.tsv"), emit: stats
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def orientationbias_options = ''
    if (orientationbias) {
        orientationbias_options = '--orientation-bias-artifact-priors ' + orientationbias.join(' --orientation-bias-artifact-priors ')
    }

    def segmentation_options = ''
    if (segmentation) {
        segmentation_options = '--tumor-segmentation ' + segmentation.join(' --tumor-segmentation ')
    }

    def contamination_options = contaminationest ? " --contamination-estimate ${contaminationest} " : ''
    if (contaminationfile) {
        contamination_options = '--contamination-table ' + contaminationfile.join(' --contamination-table ')
    }
    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK FilterMutectCalls] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" FilterMutectCalls \\
        -R $fasta \\
        -V $vcf \\
        $orientationbias_options \\
        $segmentation_options \\
        $contamination_options \\
        -O ${prefix}.vcf.gz \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
