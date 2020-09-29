process HTSLIB_TABIX {
    tag "${vcf}"

    container 'quay.io/biocontainers/tabix:0.2.6--ha92aebf_0'

    conda (params.conda ? "bioconda::tabix=0.2.6" : null)

    input:
        path vcf

    output:
        path "${vcf}.tbi"

    script:
    """
    tabix -p vcf ${vcf}
    """
}
