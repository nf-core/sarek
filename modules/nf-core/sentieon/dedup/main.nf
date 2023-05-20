process SENTIEON_DEDUP {
    tag "$meta.id"
    label 'process_medium'
    label 'sentieon'

    secret 'SENTIEON_LICENSE_BASE64'

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "Sentieon modules does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    container 'docker.io/nfcore/sentieon:202112.06'

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fasta_fai

    output:
    tuple val(meta), path("*.cram"),    emit: cram, optional: true
    tuple val(meta), path("*.crai"),    emit: crai  // Sentieon will generate a .crai AND a .bai no matter which output file type is chosen.
    tuple val(meta), path("*.bam"),     emit: bam,  optional: true
    tuple val(meta), path("*.bai"),     emit: bai
    tuple val(meta), path("*.score"),   emit: score
    tuple val(meta), path("*.metrics"), emit: metrics
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: ".cram"   // The suffix should be either ".cram" or ".bam".
    def metrics = task.ext.metrics ?: "${prefix}${suffix}.metrics"
    def sentieon_auth_mech_base64 = task.ext.sentieon_auth_mech_base64 ?: ''
    def sentieon_auth_data_base64 = task.ext.sentieon_auth_data_base64 ?: ''
    def input_list = bam.collect{"-i $it"}.join(' ')

    """
    export SENTIEON_LICENSE=\$(echo -n "\$SENTIEON_LICENSE_BASE64" | base64 -d)

    if  [ ${sentieon_auth_mech_base64} ] && [ ${sentieon_auth_data_base64} ]; then
        # If sentieon_auth_mech_base64 and sentieon_auth_data_base64 are non-empty strings, then Sentieon is mostly likely being run with some test-license.
        export SENTIEON_AUTH_MECH=\$(echo -n "${sentieon_auth_mech_base64}" | base64 -d)
        export SENTIEON_AUTH_DATA=\$(echo -n "${sentieon_auth_data_base64}" | base64 -d)
        echo "Decoded and exported Sentieon test-license system environment variables"
    fi

    sentieon driver $args $input_list -r ${fasta} --algo LocusCollector $args2 --fun score_info ${prefix}.score
    sentieon driver $args3 -t $task.cpus $input_list -r ${fasta} --algo Dedup $args4 --score_info ${prefix}.score --metrics ${metrics} ${prefix}${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cram
    touch ${prefix}.cram.crai
    touch ${prefix}.metrics
    touch ${prefix}.score

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
