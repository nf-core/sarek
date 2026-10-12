process SENTIEON_VARCAL {
    tag "${meta.id}"
    label 'process_low'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/73e9111552beb76e2ad3ad89eb75bed162d7c5b85b2433723ecb4fc96a02674a/data'
        : 'community.wave.seqera.io/library/sentieon:202503.02--def60555294d04fa'}"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path resource_vcf
    path resource_tbi
    val labels
    path fasta
    path fai

    output:
    tuple val(meta), path("*.recal"),    emit: recal
    tuple val(meta), path("*.idx"),      emit: idx
    tuple val(meta), path("*.tranches"), emit: tranches
    tuple val(meta), path("*plots.R"),   emit: plots, optional: true
    tuple val("${task.process}"), val('sentieon'), eval('sentieon driver --version | sed "s/.*-//g"'), topic: versions, emit: versions_sentieon

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Process labels to create the command string
    // labels is a list. Here is an example of what labels might look like:
    // ['--resource:dbsnp,known=false,training=true,truth=false,prior=2.0 dbsnp_146.hg38.vcf.gz', '--resource:gatk,known=false,training=true,truth=true,prior=10.0 Homo_sapiens_assembly38.known_indels.vcf.gz --resource:mills,known=false,training=true,truth=true,prior=10.0 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz']
    def labels_command = ''
    def labels_input = labels
    if (labels_input instanceof String && !labels_input.trim().isEmpty()) {
        // Process string input
        def resourceStrings = labels_input.split('--resource:').findAll()
        def processedResources = resourceStrings.collect { resource_string ->
            def items = resource_string.split(' ', 2)
            if (items.size() != 2) {
                error("Expected the resource string '${resource_string}' to contain two elements separated by a space.")
            }
            "--resource ${items[1]} --resource_param ${items[0].replaceFirst('^--resource:', '')}"
        }
        labels_command = processedResources.join(' ')
    }
    else if (labels_input instanceof List) {
        // Each list element may itself contain multiple `--resource:` directives
        // (e.g. when params.known_indels_vqsr packs both `gatk` and `mills` into
        // one string), so split each element on `--resource:` just like the
        // String branch does, then process each individual resource.
        def processedResources = labels_input.collectMany { label ->
            label.split('--resource:').findAll().collect { resource_string ->
                def items = resource_string.split(' ', 2)
                if (items.size() != 2) {
                    error("Expected the resource string '${resource_string}' to contain two elements separated by a space.")
                }
                "--resource ${items[1]} --resource_param ${items[0].replaceFirst('^--resource:', '')}"
            }
        }
        labels_command = processedResources.join(' ')
    }
    else if (labels_input != null) {
        error("Expected 'labels' to be either a String or a List, but got ${labels_input.getClass()}")
    }

    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64
        ? "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; "
        : ""
    """
    ${sentieonLicense}

    sentieon driver \\
        -r ${fasta} \\
        --algo VarCal \\
        -v ${vcf} \\
        --tranches_file ${prefix}.tranches \\
        ${labels_command} \\
        ${args} \\
        ${prefix}.recal
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.recal
    touch ${prefix}.idx
    touch ${prefix}.tranches
    touch ${prefix}plots.R
    """
}
