// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CONCAT_VCF {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::htslib=1.12" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/htslib:1.12--h9093b5e_1"
    } else {
        container "quay.io/biocontainers/htslib:1.12--h9093b5e_1"
    }

    input:
    tuple val(meta), path(vcf)
    path fai
    path bed

    output:
    tuple val(meta), path("*_*.vcf.gz"), path("*_*.vcf.gz.tbi"), emit: vcf

    script:
    def prefix           = options.suffix ? "${options.suffix}_${meta.id}" : "${meta.id}"
    def target_options   = params.target_bed ? "-t ${bed}" : ""
    def interval_options = params.no_intervals ? "-n" : ""
    """
    concatenateVCFs.sh -i ${fai} -c ${task.cpus} -o ${prefix}.vcf ${target_options} ${interval_options}
    """
}