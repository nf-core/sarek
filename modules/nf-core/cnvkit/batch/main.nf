process CNVKIT_BATCH {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-780d630a9bb6a0ff2e7b6f730906fd703e40e98f:c94363856059151a2974dc501fb07a0360cc60a3-0' :
        'biocontainers/mulled-v2-780d630a9bb6a0ff2e7b6f730906fd703e40e98f:c94363856059151a2974dc501fb07a0360cc60a3-0' }"

    input:
    tuple val(meta), path(tumor), path(normal)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    tuple val(meta4), path(targets)
    tuple val(meta5), path(reference)
    val   panel_of_normals

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

    def tumor_exists = tumor ? true : false
    def normal_exists = normal ? true : false
    def reference_exists = reference ? true : false

    // execute samtools only when cram files are input, cnvkit runs natively on bam but is prohibitively slow
    def tumor_cram = tumor_exists && tumor.Extension == "cram" ? true : false
    def normal_cram = normal_exists && normal.Extension == "cram" ? true : false
    def tumor_bam = tumor_exists && tumor.Extension == "bam" ? true : false
    def normal_bam = normal_exists && normal.Extension == "bam" ? true : false

    def tumor_out = tumor_cram ? tumor.BaseName + ".bam" : "${tumor}"

    // tumor_only mode does not need fasta & target
    // instead it requires a pre-computed reference.cnn which is built from fasta & target
    def (normal_out, normal_args, fasta_args) = ["", "", ""]
    def fai_reference = fasta_fai ? "--fai-reference ${fasta_fai}" : ""

    if (normal_exists){
        def normal_prefix = normal.BaseName
        normal_out = normal_cram ? "${normal_prefix}" + ".bam" : "${normal}"
        fasta_args = fasta ? "--fasta $fasta" : ""

        // germline mode
        // normal samples must be input without a flag
        // requires flag --normal to be empty []
        if(!tumor_exists){
            tumor_out = "${normal_prefix}" + ".bam"
            normal_args = "--normal "
        }
        // somatic mode
        else {
            normal_args = normal_prefix ? "--normal $normal_out" : ""
        }
        if (reference_exists){
            fasta_args = ""
            normal_args = ""
        }
    }

    // generation of panel of normals
    def generate_pon = panel_of_normals ? true : false

    if (generate_pon && !tumor_exists){
        def pon_input = normal.join(' ')
        normal_args = "--normal $pon_input"
        tumor_out = ""
    }

    def target_args = targets && !reference_exists ? "--targets $targets" : ""
    def reference_args = reference ? "--reference $reference" : ""

    def samtools_cram_convert = ''
    samtools_cram_convert += normal_cram ? "    samtools view -T $fasta $fai_reference $normal -@ $task.cpus -o $normal_out\n" : ''
    samtools_cram_convert += normal_cram ? "    samtools index $normal_out\n" : ''
    samtools_cram_convert += tumor_cram ? "    samtools view -T $fasta $fai_reference $tumor -@ $task.cpus -o $tumor_out\n" : ''
    samtools_cram_convert += tumor_cram ? "    samtools index $tumor_out\n" : ''
    def versions = normal_cram || tumor_cram ?
        "samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')\n        cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')" :
        "cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')"
    """
    $samtools_cram_convert

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
        ${versions}
    END_VERSIONS
    """
}
