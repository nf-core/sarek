process HTSLIB_TABIX {
    tag {vcf}

    container 'quay.io/biocontainers/tabix:0.2.6--ha92aebf_0'

    input:
        path(vcf)

    output:
        path("${vcf}.tbi")

    script:
    """
    tabix -p vcf ${vcf}
    """
}
