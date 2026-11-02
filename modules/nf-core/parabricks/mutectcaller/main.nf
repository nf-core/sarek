// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process PARABRICKS_MUTECTCALLER {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'
    // needed by the module to work properly can be removed when fixed upstream - see: https://github.com/nf-core/modules/issues/7226
    stageInMode 'copy'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.8.0-1"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(intervals)
    tuple val(meta2), path(fasta), path(fasta_fai)
    path panel_of_normals
    path panel_of_normals_index

    output:
    tuple val(meta), path("${prefix}.vcf.gz"),       emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi"),   emit: tbi
    tuple val(meta), path("${prefix}.vcf.gz.stats"), emit: stats
    path "compatible_versions.yml",                  emit: compatible_versions, optional: true
    tuple val("${task.process}"), val('parabricks'), eval("pbrun version | grep -m1 '^pbrun:' | sed 's/^pbrun:[[:space:]]*//'"), topic: versions, emit: versions_parabricks

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit(1, "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def intervals_command  = intervals     ? (intervals instanceof List ? intervals.collect { interval -> "--interval-file ${interval}" }.join(' ') : "--interval-file ${intervals}") : ""
    // prepon and postpon expect the index at <pon>.tbi
    def pon_index_command = panel_of_normals && panel_of_normals_index.name != "${panel_of_normals.name}.tbi" ? "cp -L ${panel_of_normals_index} ${panel_of_normals}.tbi" : ""
    def prepon_command = panel_of_normals ? "pbrun prepon --in-pon-file ${panel_of_normals}" : ""
    def pon_command    = panel_of_normals ? "--pon ${panel_of_normals}" : ""
    // pbrun postpon requires uncompressed .vcf input; output uncompressed when PON is provided
    def mutect_out     = panel_of_normals ? "${prefix}.vcf"    : "${prefix}.vcf.gz"
    def stats_rename   = panel_of_normals ? "mv ${prefix}.vcf.stats ${prefix}.vcf.gz.stats" : ""
    def postpon_command = panel_of_normals ? "pbrun postpon --in-vcf ${prefix}.vcf --in-pon-file ${panel_of_normals} --out-vcf ${prefix}.vcf.gz" : ""

    def num_gpus = task.accelerator ? "--num-gpus ${task.accelerator.request}" : ""
    def normal_bam_flag = normal_bam ? "--in-normal-bam ${normal_bam}" : ""
    """
    # if panel of normals specified, place its index at <pon>.tbi and run prepon
    ${pon_index_command}
    ${prepon_command}

    pbrun \\
        mutectcaller \\
        --ref ${fasta} \\
        --in-tumor-bam ${tumor_bam} \\
        ${normal_bam_flag} \\
        ${pon_command} \\
        --out-vcf ${mutect_out} \\
        ${intervals_command} \\
        ${num_gpus} \\
        ${args}

    # if panel of normals specified, rename stats and run postpon
    ${stats_rename}
    ${postpon_command}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.vcf.gz.stats

    # Capture the full version output once and store it in a variable
    pbrun_version_output=\$(pbrun mutectcaller --version 2>&1)

    # Generate compatible_versions.yml
    cat <<EOF > compatible_versions.yml
    "${task.process}":
        pbrun_version: \$(echo "\$pbrun_version_output" | grep "pbrun:" | awk '{print \$2}')
        compatible_with:
        \$(echo "\$pbrun_version_output" | awk '/Compatible With:/,/^---/{ if (\$1 ~ /^[A-Z]/ && \$1 != "Compatible" && \$1 != "---") { printf "  %s: %s\\n", \$1, \$2 } }')
    EOF
    """
}
