process CONCAT_VCF {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bcftools=1.12" : null)
    //TODO: No singularity container at the moment, use docker container for the moment
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/bcftools:1.12--h45bccc9_1' :
        'quay.io/biocontainers/bcftools:1.12--h45bccc9_1' }"

    input:
    tuple val(meta), path(vcf)
    path  fasta_fai
    path  bed

    output:
    tuple val(meta), path("*_*.vcf.gz"), path("*_*.vcf.gz.tbi"), emit: vcf

    script:
    def prefix           = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def target_options   = params.target_bed ? "-t ${bed}" : ""
    def interval_options = params.no_intervals ? "-n" : ""
    """
    concatenateVCFs.sh -i ${fasta_fai} -c ${task.cpus} -o ${prefix}.vcf ${target_options} ${interval_options}
    """
}
