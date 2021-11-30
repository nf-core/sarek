process INDEX_TARGET_BED {
    tag "$target_bed"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bcftools=1.12" : null)
    //TODO: No singularity container at the moment, use docker container for the moment
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/bcftools:1.12--h45bccc9_1' :
        'quay.io/biocontainers/bcftools:1.12--h45bccc9_1' }"

    input:
    path  target_bed

    output:
    tuple path("${target_bed}.gz"), path("${target_bed}.gz.tbi")

    script:
    """
    bgzip --threads ${task.cpus} -c ${target_bed} > ${target_bed}.gz
    tabix ${target_bed}.gz
    """
}
