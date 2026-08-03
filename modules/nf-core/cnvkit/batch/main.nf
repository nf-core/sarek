process CNVKIT_BATCH {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3e/3e8542cdb0190cfe2cedd74f714f021a2ffa94be3ec2a5b95ff52610cb3e2c34/data'
        : 'community.wave.seqera.io/library/cnvkit_htslib_samtools:86928c121163aca7'}"

    input:
    tuple val(meta), path(tumor), path(tumor_index), path(normal), path(normal_index)
    tuple val(meta2), path(fasta), path(fasta_fai)
    tuple val(meta4), path(targets)
    tuple val(meta5), path(reference)
    val panel_of_normals

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.cnn"), emit: cnn, optional: true
    tuple val(meta), path("*.cnr"), emit: cnr, optional: true
    tuple val(meta), path("*.cns"), emit: cns, optional: true
    tuple val(meta), path("*.pdf"), emit: pdf, optional: true
    tuple val(meta), path("*.png"), emit: png, optional: true
    tuple val("${task.process}"), val('cnvkit'), eval('cnvkit.py version | sed -e "s/cnvkit v//g"'), emit: versions_cnvkit, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: ''
    def tumor_exists = tumor ? true : false
    def normal_exists = normal ? true : false
    def reference_exists = reference ? true : false
    // execute samtools only when cram files are input, cnvkit runs natively on cram but is prohibitively slow
    def tumor_cram = tumor_exists && tumor.Extension == "cram" ? true : false
    def normal_cram = normal_exists && normal.Extension == "cram" ? true : false

    def tumor_out = tumor_cram ? tumor.BaseName + ".bam" : "${tumor}"

    // tumor_only mode does not need fasta & target
    // instead a pre-computed reference.cnn may be supplied which is built from fasta & target
    def (normal_out, normal_args, fasta_args) = ["", "", ""]
    def fai_reference = fasta_fai ? "--fai-reference ${fasta_fai}" : ""

    if (normal_exists) {
        def normal_prefix = normal.BaseName
        normal_out = normal_cram ? "${normal_prefix}" + ".bam" : "${normal}"
        fasta_args = fasta ? "--fasta ${fasta}" : ""

        // germline mode
        // normal samples must be input without a flag
        // requires flag --normal to be empty []
        if (!tumor_exists) {
            tumor_out = "${normal_prefix}" + ".bam"
            normal_args = "--normal "
        }
        else {
            normal_args = normal_prefix ? "--normal ${normal_out}" : ""
        }
        if (reference_exists) {
            fasta_args = ""
            normal_args = ""
        }
    }
    // generation of panel of normals
    def generate_pon = panel_of_normals ? true : false

    if (generate_pon && !tumor_exists) {
        def pon_input = normal.join(' ')
        normal_args = "--normal ${pon_input}"
        tumor_out = ""
    }

    // tumor_only mode and no reference
    // generate a "flat" reference which assumes equal coverage
    // by passing '--normal' without any files
    if (!reference_exists & !normal_exists & tumor_exists) {
        normal_args = normal_args ?: "--normal"
    }

    def target_args = targets && !reference_exists ? "--targets ${targets}" : ""
    def reference_args = reference ? "--reference ${reference}" : ""

    def samtools_cram_convert = ''
    samtools_cram_convert += normal_cram ? "    samtools view -T ${fasta} ${fai_reference} ${normal} -@ ${task.cpus} -o ${normal_out}\n" : ''
    samtools_cram_convert += normal_cram ? "    samtools index ${normal_out}\n" : ''
    samtools_cram_convert += tumor_cram ? "    samtools view -T ${fasta} ${fai_reference} ${tumor} -@ ${task.cpus} -o ${tumor_out}\n" : ''
    samtools_cram_convert += tumor_cram ? "    samtools index ${tumor_out}\n" : ''
    """
    ${samtools_cram_convert}
    cnvkit.py \\
        batch \\
        ${tumor_out} \\
        ${normal_args} \\
        ${fasta_args} \\
        ${reference_args} \\
        ${target_args} \\
        --processes ${task.cpus} \\
        ${args}
    """
    stub:
    def tumor_exists = tumor ? true : false
    def reference_exists = reference ? true : false
    // identify BED naming pattern
    def bed_prefix = reference_exists ? reference.BaseName : targets ? targets.BaseName : ""
    def bed_suffix = reference_exists ? "-tmp.bed" : ".bed"
    // execute samtools only when cram files are input, cnvkit runs natively on cram but is prohibitively slow
    def out_base_name = tumor_exists ? tumor.BaseName : normal.BaseName
    """
    touch ${bed_prefix}.antitarget${bed_suffix}
    touch ${bed_prefix}.target${bed_suffix}
    touch "reference.cnn"
    touch ${out_base_name}.antitargetcoverage.cnn
    touch ${out_base_name}.bintest.cns
    touch ${out_base_name}.call.cns
    touch ${out_base_name}.cnr
    touch ${out_base_name}.cns
    touch ${out_base_name}.targetcoverage.cnn
    """
}
