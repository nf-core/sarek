process DEEPVARIANT_RUNDEEPVARIANT {
    tag "$meta.id"
    label 'process_high'

    // FIXME Conda is not supported at the moment
    // https://github.com/bioconda/bioconda-recipes/pull/45214#issuecomment-1890937836
    // BUG https://github.com/nf-core/modules/issues/1754
    // BUG https://github.com/bioconda/bioconda-recipes/issues/30310
    container "docker.io/google/deepvariant:1.8.0"

    input:
    tuple val(meta), path(input), path(index), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(gzi)
    tuple val(meta5), path(par_bed)

    output:
    tuple val(meta), path("${prefix}.vcf.gz")      , emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi")  , emit: vcf_tbi
    tuple val(meta), path("${prefix}.g.vcf.gz")    , emit: gvcf
    tuple val(meta), path("${prefix}.g.vcf.gz.tbi"), emit: gvcf_tbi
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DEEPVARIANT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def regions = intervals ? "--regions=${intervals}" : ""
    def par_regions = par_bed ? "--par_regions_bed=${par_bed}" : ""

    """
    /opt/deepvariant/bin/run_deepvariant \\
        --ref=${fasta} \\
        --reads=${input} \\
        --output_vcf=${prefix}.vcf.gz \\
        --output_gvcf=${prefix}.g.vcf.gz \\
        ${args} \\
        ${regions} \\
        ${par_regions} \\
        --intermediate_results_dir=tmp \\
        --num_shards=${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
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
        deepvariant: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    """
}
