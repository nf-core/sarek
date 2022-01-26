process TIDDIT_SV {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::tiddit=2.12.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tiddit:2.12.2--py38h1773678_0' :
        'quay.io/biocontainers/tiddit:2.12.2--py38h1773678_0' }"

    input:
    tuple val(meta), path(input), path(index)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.vcf")        , emit: vcf
    tuple val(meta), path("*.ploidy.tab") , emit: ploidy
    tuple val(meta), path("*.signals.tab"), emit: signals
    //tuple val(meta), path("*.wig")        , emit: wig
    //tuple val(meta), path("*.gc.wig")     , emit: gc_wig

    path  "versions.yml"                  , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta != "dummy_file.txt" ? "--ref $fasta" : ""
    """
    tiddit \\
        --sv \\
        $args \\
        --bam $input \\
        $reference \\
        -o $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tiddit: \$(echo \$(tiddit 2>&1) | sed 's/^.*TIDDIT-//; s/ .*\$//')
    END_VERSIONS
    """
}
