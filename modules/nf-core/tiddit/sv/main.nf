process TIDDIT_SV {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::tiddit=3.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tiddit:3.1.0--py39h59fae87_1' :
        'quay.io/biocontainers/tiddit:3.1.0--py39h59fae87_1' }"

    input:
    tuple val(meta), path(input), path(input_index)
    path  fasta
    path  bwa_index

    output:
    tuple val(meta), path("*.vcf")         , emit: vcf
    tuple val(meta), path("*.ploidies.tab"), emit: ploidy
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bwa_command = bwa_index ? "[[ -d $bwa_index ]] && for i in $bwa_index/*; do [[ -f $fasta && ! \"\$i\" =~ .*\"$fasta\".* ]] && ln -s \$i ${fasta}.\${i##*.} || ln -s \$i .; done" : ""

    """
    $bwa_command

    tiddit \\
        --sv \\
        $args \\
        --bam $input \\
        --ref $fasta \\
        -o $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tiddit: \$(echo \$(tiddit 2>&1) | sed 's/^.*tiddit-//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf
    touch ${prefix}.ploidies.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tiddit: \$(echo \$(tiddit 2>&1) | sed 's/^.*tiddit-//; s/ .*\$//')
    END_VERSIONS
    """
}
