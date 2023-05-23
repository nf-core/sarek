process ADD_INFO_TO_VCF {
    tag "$meta.id"

    conda "anaconda::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(vcf_gz)

    output:
    tuple val(meta), path("*.added_info.vcf"), emit: vcf
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    input="input.vcf"
    output="${vcf_gz.baseName.minus(".vcf")}.added_info.vcf"
    zcat $vcf_gz > \$input
    ## Add info header lines
    grep -E "^##" \$input > \$output
    ## Add description of new INFO value
    echo '##INFO=<ID=SOURCE,Number=1,Type=String,Description="Name of vcf-file from whence the variant came">' >> \$output
    ## Add column header
    grep -E "^#CHROM" \$input >> \$output
    ## Add SOURCE value to INFO column of variant calls
    if grep -Ev "^#" \$input; then
        grep -Ev "^#" \$input | awk 'BEGIN{FS=OFS="\t"} { \$8=="." ? \$8="SOURCE=$vcf_gz" : \$8=\$8";SOURCE=$vcf_gz"; print }' >> \$output
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
