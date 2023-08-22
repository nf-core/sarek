process SENTIEON_DEDUP {
    tag "$meta.id"
    label 'process_medium'
    label 'sentieon'

    secret 'SENTIEON_LICENSE_BASE64'

    container 'nf-core/sentieon:202112.06'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)

    output:
    tuple val(meta), path("*.cram")               , emit: cram, optional: true
    tuple val(meta), path("*.crai")               , emit: crai, optional: true
    tuple val(meta), path("*.bam")                , emit: bam , optional: true
    tuple val(meta), path("*.bai")                , emit: bai
    tuple val(meta), path("*.score")              , emit: score
    tuple val(meta), path("*.metrics")            , emit: metrics
    tuple val(meta), path("*.metrics.multiqc.tsv"), emit: metrics_multiqc_tsv
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Sentieon modules do not support Conda. Please use Docker / Singularity / Podman instead."
    }
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

    sentieon driver $args $input_list -r ${fasta} --algo LocusCollector $args2 --fun score_info ${prefix}.score
    sentieon driver $args3 -t $task.cpus $input_list -r ${fasta} --algo Dedup $args4 --score_info ${prefix}.score --metrics ${metrics} ${prefix}${suffix}
    # This following tsv-file is produced in order to get a proper tsv-file with Dedup-metrics for importing in MultiQC as "custom content".
    # It should be removed once MultiQC has a module for displaying Dedup-metrics.
    head -3 ${metrics} > ${metrics}.multiqc.tsv

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
