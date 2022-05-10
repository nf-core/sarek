process CNVKIT_BATCH {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::cnvkit=0.9.9' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvkit:0.9.9--pyhdfd78af_0' :
        'quay.io/biocontainers/cnvkit:0.9.9--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(tumor), path(normal)
    path  fasta
    path  targets
    path  reference

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.cnn"), emit: cnn, optional: true
    tuple val(meta), path("*.cnr"), emit: cnr, optional: true
    tuple val(meta), path("*.cns"), emit: cns, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def normal_args = normal ? "--normal $normal" : ""
    def fasta_args = fasta ? "--fasta $fasta" : ""
    def reference_args = reference ? "--reference $reference" : ""

    def target_args = ""
    if (args.contains("--method wgs") || args.contains("-m wgs")) {
        target_args = targets ? "--targets $targets" : ""
    }
    else {
        target_args = "--targets $targets"
    }
    """
    cnvkit.py \\
        batch \\
        $tumor \\
        $normal_args \\
        $fasta_args \\
        $reference_args \\
        $target_args \\
        --processes $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$(cnvkit.py version | sed -e "s/cnvkit v//g")
    END_VERSIONS
    """
}
