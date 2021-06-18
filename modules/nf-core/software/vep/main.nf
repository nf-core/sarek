// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options = initOptions(params.options)
params.use_cache = false
use_cache = params.use_cache
params.vep_tag = "104.3.WBcel235"
vep_tag = params.vep_tag

process VEP {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::vep=104.3" : null)
    if (!use_cache) {
        container "nfcore/vep:${vep_tag}"
    } else if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/vep:104.3--pl5262h4a94de4_0"
    } else {
        container "quay.io/biocontainers/vep:104.3--pl5262h4a94de4_0"
    }

    input:
    tuple val(meta), path(vcf)
    val (vep_genome)
    val (vep_species)
    val (vep_cache_version)
    val (use_cache)
    path (vep_cache)
    val (vep_tag)

    output:
    tuple val(meta), path("*.ann.vcf"), emit: vcf
    path "*.summary.html",              emit: reports
    path "*.version.txt",               emit: version

    script:
    def software = getSoftwareName(task.process)
    prefix       = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    dir_cache    = use_cache ? "\${PWD}/${vep_cache}" : "/.vep"
    """
    mkdir ${prefix}

    vep \
        -i ${vcf} \
        -o ${prefix}.ann.vcf \
        --assembly ${vep_genome} \
        --species ${vep_species} \
        --cache \
        --cache_version ${vep_cache_version} \
        --dir_cache ${dir_cache} \
        --everything \
        --filter_common \
        --fork ${task.cpus} \
        --format vcf \
        --per_gene \
        --stats_file ${prefix}.summary.html \
        --total_length \
        --vcf

    rm -rf ${prefix}

    echo \$(vep --help 2>&1) > ${software}.version.txt
    """
}
