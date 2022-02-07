process CONCAT_VCF {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bcftools=1.14" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.14--hde04aa1_1' :
        'quay.io/biocontainers/bcftools:1.14--hde04aa1_1' }"

    input:
    tuple val(meta), path(vcf)
    path  fasta_fai
    path  target_bed

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf
    path("${prefix}.vcf.gz.tbi")             , emit: tbi
    path  "versions.yml"                     , emit: versions

    script:
    def args = task.ext.args  ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def target_options   = target_bed ? "-t ${target_bed}" : ""

    """
    concatenateVCFs.sh -i ${fasta_fai} -c ${task.cpus} -o ${prefix}.vcf ${target_options} $args
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
