// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options = initOptions(params.options)
params.use_cache = false
use_cache = params.use_cache
params.snpeff_tag = "5.0.WBcel235"
snpeff_tag = params.snpeff_tag

process SNPEFF {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::snpeff=5.0" : null)
    if (!use_cache) {
        container "nfcore/snpeff:${snpeff_tag}"
    } else if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/snpeff:5.0--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/snpeff:5.0--hdfd78af_1"
    }

    input:
    tuple val(meta), path(vcf)
    val (snpeff_db)
    val (use_cache)
    path (snpeff_cache)
    val (snpeff_tag)

    output:
    tuple val(meta), path("*.ann.vcf"), emit: vcf
    path "*.csv",                       emit: reports
    path "*.version.txt",               emit: version

    script:
    def software = getSoftwareName(task.process)
    prefix       = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    cache        = use_cache ? "-dataDir \${PWD}/${snpeff_cache}" : ""
    """
    snpEff -Xmx${task.memory.toGiga()}g \
        ${snpeff_db} \
        -csvStats ${prefix}.csv \
        ${cache} \
        -canon \
        -v \
        ${vcf} \
        > ${prefix}.ann.vcf

    echo \$(snpEff -version 2>&1) > ${software}.version.txt
    """
}
