process CONCAT_VCF {
    label 'cpus_8'

    tag "${options.publish_dir}-${meta.id}"

    publishDir "${params.outdir}/VariantCalling/${meta.id}/${options.publish_dir}", mode: params.publish_dir_mode

    container "quay.io/biocontainers/htslib:1.11--hd3b49d5_0"

    conda (params.conda ? "bioconda::htslib=1.11" : null)

    input:
        tuple val(meta), path(vcf)
        path fai
        path bed
        val options

    output:
        tuple val(meta), path("*_*.vcf.gz"), path("*_*.vcf.gz.tbi"), emit: vcf

    script:
    name = options.suffix ? "${options.publish_dir}_${meta.id}${options.suffix}" : "${options.publish_dir}_${meta.id}"
    target_options = params.target_bed ? "-t ${bed}" : ""
    interval_options = params.no_intervals ? "-n" : ""
    """
    concatenateVCFs.sh -i ${fai} -c ${task.cpus} -o ${name}.vcf ${target_options} ${interval_options}
    """
}