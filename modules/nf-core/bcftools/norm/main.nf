process BCFTOOLS_NORM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5acacb55c52bec97c61fd34ffa8721fce82ce823005793592e2a80bf71632cd0/data':
        'community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.{vcf,vcf.gz,bcf,bcf.gz}"), emit: vcf
    tuple val(meta), path("*.tbi")                    , emit: tbi, optional: true
    tuple val(meta), path("*.csi")                    , emit: csi, optional: true
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--output-type z'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--output-type b") || args.contains("-Ob") ? "bcf.gz" :
                    args.contains("--output-type u") || args.contains("-Ou") ? "bcf" :
                    args.contains("--output-type z") || args.contains("-Oz") ? "vcf.gz" :
                    args.contains("--output-type v") || args.contains("-Ov") ? "vcf" :
                    "vcf.gz"

    """
    bcftools norm \\
        --fasta-ref ${fasta} \\
        --output ${prefix}.${extension} \\
        $args \\
        --threads $task.cpus \\
        ${vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: '--output-type z'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--output-type b") || args.contains("-Ob") ? "bcf.gz" :
                    args.contains("--output-type u") || args.contains("-Ou") ? "bcf" :
                    args.contains("--output-type z") || args.contains("-Oz") ? "vcf.gz" :
                    args.contains("--output-type v") || args.contains("-Ov") ? "vcf" :
                    "vcf.gz"
    def index = ''
    if (extension in ['vcf.gz', 'bcf', 'bcf.gz']) {
        if (['--write-index=tbi', '-W=tbi'].any { args.contains(it) }  && extension == 'vcf.gz') {
            index = 'tbi'
        } else if (['--write-index=tbi', '-W=tbi', '--write-index=csi', '-W=csi', '--write-index', '-W'].any { args.contains(it) }) {
            index = 'csi'
        }
    }
    def create_cmd = extension.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    def create_index = index ? "touch ${prefix}.${extension}.${index}" : ""

    """
    ${create_cmd} ${prefix}.${extension}
    ${create_index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
