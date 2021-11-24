// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MANTA_SOMATIC {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::manta=1.6.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/manta:1.6.0--h9ee0642_1"
    } else {
        container "quay.io/biocontainers/manta:1.6.0--h9ee0642_1"
    }

    input:
    tuple val(meta), path(input_normal), path(input_index_normal), path(input_tumor), path(input_index_tumor)
    path fasta
    path fai
    path target_bed
    path target_bed_tbi

    output:
    tuple val(meta), path("*.candidate_small_indels.vcf.gz")     , emit: candidate_small_indels_vcf
    tuple val(meta), path("*.candidate_small_indels.vcf.gz.tbi") , emit: candidate_small_indels_vcf_tbi
    tuple val(meta), path("*.candidate_sv.vcf.gz")               , emit: candidate_sv_vcf
    tuple val(meta), path("*.candidate_sv.vcf.gz.tbi")           , emit: candidate_sv_vcf_tbi
    tuple val(meta), path("*.diploid_sv.vcf.gz")                 , emit: diploid_sv_vcf
    tuple val(meta), path("*.diploid_sv.vcf.gz.tbi")             , emit: diploid_sv_vcf_tbi
    tuple val(meta), path("*.somatic_sv.vcf.gz")                 , emit: somatic_sv_vcf
    tuple val(meta), path("*.somatic_sv.vcf.gz.tbi")             , emit: somatic_sv_vcf_tbi
    path "versions.yml"                                          , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def options_manta = target_bed ? "--exome --callRegions $target_bed" : ""

    """
    configManta.py \
        --tumorBam $input_tumor \
        --normalBam $input_normal \
        --reference $fasta \
        $options_manta \
        --runDir manta

    python manta/runWorkflow.py -m local -j $task.cpus

    mv manta/results/variants/candidateSmallIndels.vcf.gz     ${prefix}.candidate_small_indels.vcf.gz
    mv manta/results/variants/candidateSmallIndels.vcf.gz.tbi ${prefix}.candidate_small_indels.vcf.gz.tbi
    mv manta/results/variants/candidateSV.vcf.gz              ${prefix}.candidate_sv.vcf.gz
    mv manta/results/variants/candidateSV.vcf.gz.tbi          ${prefix}.candidate_sv.vcf.gz.tbi
    mv manta/results/variants/diploidSV.vcf.gz                ${prefix}.diploid_sv.vcf.gz
    mv manta/results/variants/diploidSV.vcf.gz.tbi            ${prefix}.diploid_sv.vcf.gz.tbi
    mv manta/results/variants/somaticSV.vcf.gz                ${prefix}.somatic_sv.vcf.gz
    mv manta/results/variants/somaticSV.vcf.gz.tbi            ${prefix}.somatic_sv.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( configManta.py --version )
    END_VERSIONS
    """
}
