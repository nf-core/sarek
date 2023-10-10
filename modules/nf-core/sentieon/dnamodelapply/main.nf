process SENTIEON_DNAMODELAPPLY {
    tag "$meta.id"
    label 'process_high'
    label 'sentieon'

    secret 'SENTIEON_LICENSE_BASE64'

    container 'nf-core/sentieon:202112.06'

    input:
    tuple val(meta), path(vcf), path(idx)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(ml_model)

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: index
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Sentieon modules do not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args      = task.ext.args   ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def sentieon_auth_mech_base64 = task.ext.sentieon_auth_mech_base64 ?: ''
    def sentieon_auth_data_base64 = task.ext.sentieon_auth_data_base64 ?: ''

    """
    if [ "\${#SENTIEON_LICENSE_BASE64}" -lt "1500" ]; then  # If the string SENTIEON_LICENSE_BASE64 is short, then it is an encrypted url.
        export SENTIEON_LICENSE=\$(echo -e "\$SENTIEON_LICENSE_BASE64" | base64 -d)
    else  # Localhost license file
        # The license file is stored as a nextflow variable like, for instance, this:
        # nextflow secrets set SENTIEON_LICENSE_BASE64 \$(cat <sentieon_license_file.lic> | base64 -w 0)
        export SENTIEON_LICENSE=\$(mktemp)
        echo -e "\$SENTIEON_LICENSE_BASE64" | base64 -d > \$SENTIEON_LICENSE
    fi

    if  [ ${sentieon_auth_mech_base64} ] && [ ${sentieon_auth_data_base64} ]; then
        # If sentieon_auth_mech_base64 and sentieon_auth_data_base64 are non-empty strings, then Sentieon is mostly likely being run with some test-license.
        export SENTIEON_AUTH_MECH=\$(echo -n "${sentieon_auth_mech_base64}" | base64 -d)
        export SENTIEON_AUTH_DATA=\$(echo -n "${sentieon_auth_data_base64}" | base64 -d)
        echo "Decoded and exported Sentieon test-license system environment variables"
    fi

    sentieon driver \\
        -t $task.cpus \\
        -r $fasta \\
        $args \\
        --algo DNAModelApply \\
        --model $ml_model \\
        -v $vcf \\
        ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Sentieon modules do not support Conda. Please use Docker / Singularity / Podman instead."
    }
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g" )
    END_VERSIONS
    """
}
