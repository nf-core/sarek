process PARABRICKS_DEEPSOMATIC {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'
    // needed by the module to work properly can be removed when fixed upstream - see: https://github.com/nf-core/modules/issues/7226
    stageInMode 'copy'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.7.1-1"

    input:
    tuple val(meta), path(input_tumor), path(index_tumor), path(input_normal), path(index_normal), path(intervals)
    tuple val(ref_meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf, optional: true
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: vcf_tbi, optional: true
    tuple val(meta), path("${prefix}.g.vcf.gz"), emit: gvcf, optional: true
    tuple val(meta), path("${prefix}.g.vcf.gz.tbi"), emit: gvcf_tbi, optional: true
    path "compatible_versions.yml", emit: compatible_versions, optional: true
    tuple val("${task.process}"), val('parabricks'), eval("pbrun version | grep -m1 '^pbrun:' | sed 's/^pbrun:[[:space:]]*//'"), topic: versions, emit: versions_parabricks

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit(1, "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args: String = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def output_mode: String = task.ext.output_mode ?: 'vcf'
    if (!(output_mode in ['vcf', 'gvcf'])) {
        error("PARABRICKS_DEEPSOMATIC output mode must be 'vcf' or 'gvcf'.")
    }
    def output_file: String = output_mode == 'gvcf' ? "${prefix}.g.vcf.gz" : "${prefix}.vcf.gz"
    def interval_command: String = intervals ? intervals.collect { interval -> "--interval-file ${interval}" }.join(' ') : ""
    def num_gpus: String = task.accelerator ? "--num-gpus ${task.accelerator.request}" : ''
    """
    pbrun \\
        deepsomatic \\
        --ref ${fasta} \\
        --in-tumor-bam ${input_tumor} \\
        --in-normal-bam ${input_normal} \\
        --out-variants ${output_file} \\
        ${interval_command} \\
        ${num_gpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def output_mode: String = task.ext.output_mode ?: 'vcf'
    if (!(output_mode in ['vcf', 'gvcf'])) {
        error("PARABRICKS_DEEPSOMATIC output mode must be 'vcf' or 'gvcf'.")
    }
    def output_file: String = output_mode == 'gvcf' ? "${prefix}.g.vcf.gz" : "${prefix}.vcf.gz"
    """
    echo '' | gzip > ${output_file}
    touch ${output_file}.tbi

    # Capture the full version output once and store it in a variable
    pbrun_version_output=\$(pbrun deepsomatic --version 2>&1)

    # Generate compatible_versions.yml
    cat <<EOF > compatible_versions.yml
    "${task.process}":
        pbrun_version: \$(echo "\$pbrun_version_output" | grep "pbrun:" | awk '{print \$2}')
        compatible_with:
        \$(echo "\$pbrun_version_output" | awk '/Compatible With:/,/^---/{ if (\$1 ~ /^[A-Z]/ && \$1 != "Compatible" && \$1 != "---") { printf "  %s: %s\\n", \$1, \$2 } }')
    EOF
    """
}
