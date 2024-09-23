process DEEPVARIANT_POSTPROCESSVARIANTS {
    tag "$meta.id"
    label 'process_medium'

    //Conda is not supported at the moment
    container "nf-core/deepvariant:1.6.1"

    input:
    tuple val(meta), path(variant_calls_tfrecord_files), path(gvcf_tfrecords)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(gzi)

    output:
    tuple val(meta), path("${prefix}.vcf.gz")      ,  emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi")  ,  emit: vcf_tbi
    tuple val(meta), path("${prefix}.g.vcf.gz")    ,  emit: gvcf
    tuple val(meta), path("${prefix}.g.vcf.gz.tbi"),  emit: gvcf_tbi

    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DEEPVARIANT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def variant_calls_tfrecord_name = variant_calls_tfrecord_files[0].name.replaceFirst(/-\d{5}-of-\d{5}/, "")

    def gvcf_matcher = gvcf_tfrecords[0].baseName =~ /^(.+)-\d{5}-of-(\d{5})$/
    if (!gvcf_matcher.matches()) {
        throw new IllegalArgumentException("tfrecord baseName '" + gvcf_tfrecords[0].baseName + "' doesn't match the expected pattern")
    }
    def gvcf_tfrecord_name = gvcf_matcher[0][1]
    def gvcf_shardCount = gvcf_matcher[0][2]
    // Reconstruct the logical name - ${tfrecord_name}.examples.tfrecord@${task.cpus}.gz
    def gvcf_tfrecords_logical_name = "${gvcf_tfrecord_name}@${gvcf_shardCount}.gz"

    """
    /opt/deepvariant/bin/postprocess_variants \\
        ${args} \\
        --ref "${fasta}" \\
        --infile "${variant_calls_tfrecord_name}" \\
        --outfile "${prefix}.vcf.gz" \\
        --nonvariant_site_tfrecord_path "${gvcf_tfrecords_logical_name}" \\
        --gvcf_outfile "${prefix}.g.vcf.gz" \\
        --cpus $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant_postprocessvariants: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DEEPVARIANT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    echo "" | gzip > ${prefix}.g.vcf.gz
    touch ${prefix}.g.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant_callvariants: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    """
}
