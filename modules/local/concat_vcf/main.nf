process CONCAT_VCF {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bcftools=1.12" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        //TODO: No singularity container at the moment, use docker container for the moment
        container "quay.io/biocontainers/bcftools:1.12--h45bccc9_1"
    } else {
        container "quay.io/biocontainers/bcftools:1.12--h45bccc9_1"
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
