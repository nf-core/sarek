process GATK4_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gatk4=4.2.3.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.3.0--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.3.0--hdfd78af_0"
    }

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*.bam"), path("*.bai"), emit: bam_bai
    tuple val(meta), path("*.metrics")           , emit: metrics
    path "versions.yml"                          , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def bam_list = bams.collect(){ bam -> "--INPUT ".concat(bam.toString()) }.join(" ")
    def avail_mem       = 3
    if (!task.memory) {
        log.info '[GATK HaplotypeCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk MarkDuplicates \\
        $bam_list \\
        --METRICS_FILE ${prefix}.metrics \\
        --TMP_DIR . \\
        --CREATE_INDEX true \\
        --OUTPUT ${prefix}.bam \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
