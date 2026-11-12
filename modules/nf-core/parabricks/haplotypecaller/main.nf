// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process PARABRICKS_HAPLOTYPECALLER {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.8.0-1"

    input:
    tuple val(meta), path(input), path(input_index), path(intervals)
    tuple val(ref_meta), path(fasta)
    val out_ext

    output:
    tuple val(meta), path("${prefix}.vcf"),          emit: vcf,                 optional: true
    tuple val(meta), path("${prefix}.vcf.gz"),       emit: vcf_gz,              optional: true
    tuple val(meta), path("${prefix}.vcf.gz.tbi"),   emit: vcf_tbi,             optional: true
    tuple val(meta), path("${prefix}.g.vcf.gz"),     emit: gvcf,                optional: true
    tuple val(meta), path("${prefix}.g.vcf.gz.tbi"), emit: gvcf_tbi,            optional: true
    path "compatible_versions.yml",                  emit: compatible_versions, optional: true
    tuple val("${task.process}"), val('parabricks'), eval("pbrun version | grep -m1 '^pbrun:' | sed 's/^pbrun:[[:space:]]*//'"), topic: versions, emit: versions_parabricks

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit(1, "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args              = task.ext.args ?: ''
    prefix                = task.ext.prefix ?: "${meta.id}"
    def gvcf              = out_ext == 'g.vcf.gz' ? "--gvcf" : ''
    def intervals_command = intervals ? (intervals instanceof List ? intervals.collect { interval -> "--interval-file ${interval}" }.join(' ') : "--interval-file ${intervals}") : ""
    def num_gpus          = task.accelerator ? "--num-gpus ${task.accelerator.request}" : ''
    if (out_ext !in ['vcf', 'vcf.gz', 'g.vcf.gz']) {
        error("parabricks/haplotypecaller: Invalid out_ext: ${out_ext}. Allowed values are: vcf, vcf.gz, g.vcf.gz")
    }
    if (args.tokenize().contains('--gvcf')) {
        error("parabricks/haplotypecaller: Do not pass --gvcf in ext.args, set out_ext to 'g.vcf.gz' instead")
    }
    """
    pbrun \\
        haplotypecaller \\
        --ref ${fasta} \\
        --in-bam ${input} \\
        --out-variants ${prefix}.${out_ext} \\
        ${gvcf} \\
        ${intervals_command} \\
        ${num_gpus} \\
        --preserve-file-symlinks \\
        ${args}
    """

    stub:
    prefix         = task.ext.prefix ?: "${meta.id}"
    def output_cmd = out_ext.endsWith(".gz") ? "echo \"\" | gzip > ${prefix}.${out_ext}; touch ${prefix}.${out_ext}.tbi" : "touch ${prefix}.${out_ext}"
    """
    ${output_cmd}

    # Capture the full version output once and store it in a variable
    pbrun_version_output=\$(pbrun haplotypecaller --version 2>&1)

    # Generate compatible_versions.yml
    cat <<EOF > compatible_versions.yml
    "${task.process}":
        pbrun_version: \$(echo "\$pbrun_version_output" | grep "pbrun:" | awk '{print \$2}')
        compatible_with:
        \$(echo "\$pbrun_version_output" | awk '/Compatible With:/,/^---/{ if (\$1 ~ /^[A-Z]/ && \$1 != "Compatible" && \$1 != "---") { printf "  %s: %s\\n", \$1, \$2 } }')
    EOF
    """
}
