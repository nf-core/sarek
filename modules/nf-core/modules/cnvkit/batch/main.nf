process CNVKIT_BATCH {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::cnvkit=0.9.9 bioconda::samtools=1.15.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-780d630a9bb6a0ff2e7b6f730906fd703e40e98f:304d1c5ab610f216e77c61420ebe85f1e7c5968a-0' :
        'quay.io/biocontainers/mulled-v2-780d630a9bb6a0ff2e7b6f730906fd703e40e98f:304d1c5ab610f216e77c61420ebe85f1e7c5968a-0' }"

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
    tuple val(meta), path("*.pdf"), emit: pdf, optional: true
    tuple val(meta), path("*.png"), emit: png, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    // execute samtools only when cram files are input, cnvkit runs natively on bam but is prohibitively slow
    // input pair is assumed to have same extension if both exist
    def is_cram = tumor.Extension == "cram" ? true : false
    def tumor_out = is_cram ? tumor.BaseName + ".bam" : "${tumor}"

    // do not run samtools on normal samples in tumor_only mode
    def normal_exists = normal ? true: false
    // tumor_only mode does not need fasta & target
    // instead it requires a pre-computed reference.cnn which is built from fasta & target
    def (normal_out, normal_args, fasta_args) = ["", "", ""]

    if (normal_exists){
        def normal_prefix = normal.BaseName
        normal_out = is_cram ? "${normal_prefix}" + ".bam" : "${normal}"
        normal_args = normal_prefix ? "--normal $normal_out" : ""
        fasta_args = fasta ? "--fasta $fasta" : ""
    }

    def target_args = targets ? "--targets $targets" : ""
    def reference_args = reference ? "--reference $reference" : ""

    """
    if $is_cram; then
        samtools view -T $fasta $tumor -@ $task.cpus -o $tumor_out
        if $normal_exists; then
            samtools view -T $fasta $normal -@ $task.cpus -o $normal_out
        fi
    fi

    cnvkit.py \\
        batch \\
        $tumor_out \\
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
