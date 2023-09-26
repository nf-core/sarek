process SENTIEON_VARCAL {
    tag "$meta.id"
    label 'process_low'
    label 'sentieon'

    secret 'SENTIEON_LICENSE_BASE64'

    container 'nf-core/sentieon:202112.06'

    input:
    tuple val(meta), path(vcf), path(tbi) // input vcf and tbi of variants to recalibrate
    path resource_vcf   // resource vcf
    path resource_tbi   // resource tbi
    val labels          // string (or list of strings) containing dedicated resource labels already formatted with '--resource:' tag
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.recal")   , emit: recal
    tuple val(meta), path("*.idx")     , emit: idx
    tuple val(meta), path("*.tranches"), emit: tranches
    tuple val(meta), path("*plots.R")  , emit: plots, optional:true
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Sentieon modules do not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference_command = fasta ? "--reference $fasta " : ''
    def labels_command = ''

    // labels is a list. Here is an example of what labels might look like:
    // ['--resource:dbsnp,known=false,training=true,truth=false,prior=2.0 dbsnp_146.hg38.vcf.gz', '--resource:gatk,known=false,training=true,truth=true,prior=10.0 Homo_sapiens_assembly38.known_indels.vcf.gz --resource:mills,known=false,training=true,truth=true,prior=10.0 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz']
    for(label in labels){
        for(gatk_resource_string in label.split('--resource:').findAll()){  // The findAll cmd is there to remove any empty string elements in the list
            def items = gatk_resource_string.split(' ')
            // Here is an example of what the list items might look like:
            // ['dbsnp,known=false,training=true,truth=false,prior=2.0', 'dbsnp_146.hg38.vcf.gz']
            if (items.size() != 2) {
                error("Expected the list '${items}' to contain two elements.")
            }
            labels_command +=  "--resource ${items[1]} --resource_param ${items[0]} "
        }
    }

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

    sentieon driver -r ${fasta}  --algo VarCal \\
        -v $vcf \\
        --tranches_file ${prefix}.tranches \\
        $labels_command \\
        $args \\
        ${prefix}.recal

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Sentieon modules do not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix   = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.recal
    touch ${prefix}.idx
    touch ${prefix}.tranches
    touch ${prefix}plots.R

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
