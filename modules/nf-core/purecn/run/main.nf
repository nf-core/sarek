process PURECN_RUN {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-582ac26068889091d5e798347c637f8208d77a71:a29c64a63498b1ee8b192521fdf6ed3c65506994-0':
        'biocontainers/mulled-v2-582ac26068889091d5e798347c637f8208d77a71:a29c64a63498b1ee8b192521fdf6ed3c65506994-0' }"

    input:
    tuple val(meta), path(intervals), path(coverage)
    path normal_db
    val genome

    output:
    tuple val(meta), path("*.pdf")                             , emit: pdf
    tuple val(meta), path("*_local_optima.pdf")                , emit: local_optima_pdf
    tuple val(meta), path("*_dnacopy.seg")                     , emit: seg
    tuple val(meta), path("*_genes.csv")                       , emit: genes_csv                   , optional: true
    tuple val(meta), path("*_amplification_pvalues.csv")       , emit: amplification_pvalues_csv   , optional: true
    tuple val(meta), path("*.vcf.gz")                          , emit: vcf_gz                      , optional: true
    tuple val(meta), path("*_variants.csv")                    , emit: variants_csv                , optional: true
    tuple val(meta), path("*_loh.csv")                         , emit: loh_csv                     , optional: true
    tuple val(meta), path("*_chromosomes.pdf")                 , emit: chr_pdf                     , optional: true
    tuple val(meta), path("*_segmentation.pdf")                , emit: segmentation_pdf            , optional: true
    tuple val(meta), path("*_multisample.seg")                 , emit: multisample_seg             , optional: true
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.4.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    library_path=\$(Rscript -e 'cat(.libPaths(), sep = "\\n")')
    Rscript "\$library_path"/PureCN/extdata/PureCN.R \\
        --out ./ \\
        --tumor ${coverage} \\
        --sampleid ${prefix} \\
        --normaldb ${normal_db} \\
        --intervals ${intervals} \\
        --genome ${genome} \\
        --parallel \\
        --cores ${task.cpus} \\
        --stats-file ${prefix}_stats.txt \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purecn: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.4.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch ${prefix}.pdf
    touch ${prefix}_local_optima.pdf
    touch ${prefix}_dnacopy.seg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purecn: ${VERSION}
    END_VERSIONS
    """
}
