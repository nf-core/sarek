process MANTA_SOMATIC {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::manta=1.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/manta:1.6.0--h9ee0642_1' :
        'quay.io/biocontainers/manta:1.6.0--h9ee0642_1' }"

    input:
    tuple val(meta), path(input_normal), path(input_index_normal), path(input_tumor), path(input_index_tumor), path(target_bed), path(target_bed_tbi)
    path fasta
    path fai

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

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
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
    "${task.process}":
        manta: \$( configManta.py --version )
    END_VERSIONS
    """
}
