// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SAMTOOLS_MERGE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::samtools=1.14" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0"
    } else {
        container "quay.io/biocontainers/samtools:1.14--hb421002_0"
    }

    input:
    tuple val(meta), path(input_files)
    path fasta

    output:
    tuple val(meta), path("${prefix}.bam"),  optional:true, emit: bam
    tuple val(meta), path("${prefix}.cram"), optional:true, emit: cram
    path  "versions.yml"                                  , emit: versions

    script:
    prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def file_type = input_files[0].getExtension()
    def reference = fasta ? "--reference ${fasta}" : ""
    """
    samtools merge --threads ${task.cpus-1} $options.args ${reference} ${prefix}.${file_type} $input_files

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
