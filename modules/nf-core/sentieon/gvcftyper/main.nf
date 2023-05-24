process SENTIEON_GVCFTYPER {
    tag "$meta.id"
    label 'process_high'
    label 'sentieon'

    secret 'SENTIEON_LICENSE_BASE64'

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "Sentieon modules does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    container 'docker.io/nfcore/sentieon:202112.06'

    input:
    tuple val(meta), path(gvcfs), path(tbis), path(intervals)
    path  fasta
    path  fai
    path  dbsnp
    path  dbsnp_tbi

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf_gz
    tuple val(meta), path("*.vcf.gz.tbi"), emit: vcf_gz_tbi
    path("versions.yml")                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sentieon_auth_mech_base64 = task.ext.sentieon_auth_mech_base64 ?: ''
    def sentieon_auth_data_base64 = task.ext.sentieon_auth_data_base64 ?: ''
    def gvcfs_input = '-v ' + gvcfs.join(' -v ')
    def dbsnp_cmd = dbsnp ? "--dbsnp $dbsnp" : ""

    """
    export SENTIEON_LICENSE=\$(echo -n "\$SENTIEON_LICENSE_BASE64" | base64 -d)

    if  [ ${sentieon_auth_mech_base64} ] && [ ${sentieon_auth_data_base64} ]; then
        # If sentieon_auth_mech_base64 and sentieon_auth_data_base64 are non-empty strings, then Sentieon is mostly likely being run with some test-license.
        export SENTIEON_AUTH_MECH=\$(echo -n "${sentieon_auth_mech_base64}" | base64 -d)
        export SENTIEON_AUTH_DATA=\$(echo -n "${sentieon_auth_data_base64}" | base64 -d)
        echo "Decoded and exported Sentieon test-license system environment variables"
    fi

    sentieon driver -r ${fasta} --algo GVCFtyper ${gvcfs_input} ${dbsnp_cmd} ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
