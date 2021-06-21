// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_APPLYBQSR_SPARK {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gatk4==4.1.9.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.1.9.0--py39_0"
    } else {
        container "quay.io/biocontainers/gatk4:4.1.9.0--py39_0"
    }

    input:
    tuple val(meta), path(cram), path(crai), path(bqsr_table), path(intervalsBed)
    path fasta
    path fastaidx
    path dict

    output:
    tuple val(meta), path("*.cram"), emit: cram
    path "*.version.txt",            emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def intervalsCommand = intervalsBed ? "-L ${intervalsBed}" : ""
    """
    gatk ApplyBQSRSpark \\
        -R $fasta \\
        -I $cram \\
        --bqsr-recal-file $bqsr_table \\
        $intervalsCommand \\
        --tmp-dir . \
        -O ${prefix}.cram \\
        $options.args \
        --spark-master local[${task.cpus}]


    gatk --version | grep Picard | sed "s/Picard Version: //g" > ${software}.version.txt
    """
}
