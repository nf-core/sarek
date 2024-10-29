process PURECN_INTERVALFILE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-org.hs.eg.db_bioconductor-purecn_bioconductor-txdb.hsapiens.ucsc.hg19.knowngene_bioconductor-txdb.hsapiens.ucsc.hg38.knowngene_pruned:c04754ed02eb7cd3':
        'community.wave.seqera.io/library/bioconductor-org.hs.eg.db_bioconductor-purecn_bioconductor-txdb.hsapiens.ucsc.hg19.knowngene_bioconductor-txdb.hsapiens.ucsc.hg38.knowngene_pruned:cd51f6d3c90eb24f' }"

    input:
    tuple val(meta), path(target_bed)
    tuple val(meta2), path(fasta)
    val   genome

    output:
    tuple val(meta), path("*.txt"), emit: txt
    // Only produced if --export is used
    tuple val(meta), path("*.bed"), emit: bed, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    library_path=\$(Rscript -e 'cat(.libPaths(), sep = "\\n")')
    Rscript "\$library_path"/PureCN/extdata/IntervalFile.R \\
        --in-file ${target_bed} \\
        --fasta ${fasta} \\
        --out-file ${prefix}.txt \\
        --genome ${genome} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purecn: \$(Rscript -e 'packageVersion("PureCN")' | sed -n 's|\\[1\\] ‘\\(.*\\)’|\\1|p')
    END_VERSIONS
    """

    stub:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bed = args.contains("--export") ? "touch ${prefix}.bed" : ""

    """
    touch ${prefix}.txt
    ${bed}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purecn: \$(Rscript -e 'packageVersion("PureCN")' | sed -n 's|\\[1\\] ‘\\(.*\\)’|\\1|p')
    END_VERSIONS
    """
}
