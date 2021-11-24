// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_APPLYBQSR {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.0.0--0"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.0.0--0"
    }

    input:
    tuple val(meta), path(cram), path(crai), path(bqsr_table), path(intervalsBed)
    path(fasta)
    path(fastaidx)
    path(dict)

    output:
    tuple val(meta), path("*.cram"), emit: cram
    path "versions.yml"            , emit: versions

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def intervalsCommand = intervalsBed ? "-L ${intervalsBed}" : ""
    if (!task.memory) {
        log.info '[GATK ApplyBQSR] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk ApplyBQSR \\
        -R $fasta \\
        -I $cram \\
        --bqsr-recal-file $bqsr_table \\
        $intervalsCommand \\
        --tmp-dir . \\
        -O ${prefix}.cram \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
