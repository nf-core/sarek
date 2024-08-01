process SENTIEON_BWAMEM {
    tag "$meta.id"
    label 'process_high'
    label 'sentieon'

    // NOTE this is not required for the process to run, but it is required for the process to be run in GitHub actions or nf-core MegaTests
    // The rest of the secrets aren't really "secrets" because they're not sensitive information for users
    // secret SENTIEON_AUTH_DATA

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sentieon:202308.02--h43eeafb_0' :
        'biocontainers/sentieon:202308.02--h43eeafb_0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fasta_fai)

    output:
    tuple val(meta), path("*.bam"), path("*.bai"), emit: bam_and_bai
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The following code sets LD_LIBRARY_PATH in the script-section when the module is run by Singularity.
    // That turned out to be one way of overcoming the following issue with the Singularity-Sentieon-containers from galaxy, Sentieon (LD_LIBRARY_PATH) and the way Nextflow runs Singularity-containers.
    // The galaxy container uses a runscript which is responsible for setting LD_PRELOAD properly. Nextflow executes singularity containers using `singularity exec`, which avoids the run script, leading to the LD_LIBRARY_PATH/libstdc++.so.6 error.
    if (workflow.containerEngine in ['singularity','apptainer']) {
        fix_ld_library_path = 'LD_LIBRARY_PATH=/usr/local/lib/:\$LD_LIBRARY_PATH;export LD_LIBRARY_PATH'
    } else {
        fix_ld_library_path = ''
    }

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sentieon_auth_mech_base64 = task.ext.sentieon_auth_mech_base64 ?: ''
    def sentieon_auth_data_base64 = task.ext.sentieon_auth_data_base64 ?: ''

    """
    if [ "\${SENTIEON_LICSRVR_IP}" ]; then
        # NOTE: This is how pipeline users will use Sentieon in real world
        echo "Using a Sentieon License Server"
        export SENTIEON_LICENSE="\${SENTIEON_LICSRVR_IP}"
    else
        # NOTE: This is how pipeline users will test out Sentieon
        echo "Localhost license file"
        # The license file is stored as a nextflow variable like, for instance, this:
        # nextflow secrets set SENTIEON_LICENSE_BASE64 \$(cat <sentieon_license_file.lic> | base64 -w 0)
        export SENTIEON_LICENSE=\$(mktemp)
        echo -e "\$SENTIEON_LICENSE_BASE64" | base64 -d > \$SENTIEON_LICENSE
    fi

    if  [ "\${SENTIEON_AUTH_MECH}" ] && [ "\${SENTIEON_AUTH_DATA}" ]; then
        # NOTE: This should only happen in GitHub Actions or nf-core/megatests
        echo "If sentieon_auth_mech and sentieon_auth_data are non-empty strings, then Sentieon is mostly likely being run with some test-license."
    fi

    $fix_ld_library_path

    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

    sentieon bwa mem \\
        $args \\
        -t $task.cpus \\
        \$INDEX \\
        $reads \\
        | sentieon util sort -r $fasta -t $task.cpus -o ${prefix}.bam --sam2bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """

    stub:
    // The following code sets LD_LIBRARY_PATH in the script-section when the module is run by Singularity.
    // That turned out to be one way of overcoming the following issue with the Singularity-Sentieon-containers from galaxy, Sentieon (LD_LIBRARY_PATH) and the way Nextflow runs Singularity-containers.
    // The galaxy container uses a runscript which is responsible for setting LD_PRELOAD properly. Nextflow executes singularity containers using `singularity exec`, which avoids the run script, leading to the LD_LIBRARY_PATH/libstdc++.so.6 error.
    if (workflow.containerEngine in ['singularity','apptainer']) {
        fix_ld_library_path = 'LD_LIBRARY_PATH=/usr/local/lib/:\$LD_LIBRARY_PATH;export LD_LIBRARY_PATH'
    } else {
        fix_ld_library_path = ''
    }

    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    $fix_ld_library_path

    touch ${prefix}.bam
    touch ${prefix}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}
