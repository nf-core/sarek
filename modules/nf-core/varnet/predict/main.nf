process VARNET_PREDICT {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/varnet:1.5.3--pyhdfd78af_0':
        'quay.io/biocontainers/varnet:1.5.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(input_tumor), path(index_tumor), path(input_normal), path(index_normal), path(candidates, stageAs: 'candidates')
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("${prefix}/${prefix}.vcf.gz"), emit: vcf
    tuple val(meta), path("${prefix}/predictions"), emit: predictions
    tuple val("${task.process}"), val("varnet"), val("1.5.3"), emit: versions_varnet, topic: versions
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!input_normal) {
        error "VARNET_PREDICT requires a matched normal BAM. Tumor-only mode needs large germline resource files (dbSNP and gnomAD) that are not shipped with the conda package, so it is not supported by this module. To run tumor-only, use the VarNet Docker image directly: https://github.com/skandlab/VarNet"
    }
    """
    export TF_CPP_MIN_LOG_LEVEL=3

    # varnet-predict expects the candidates from varnet-filter at <output_dir>/<sample_name>/candidates
    mkdir -p ${prefix}
    cp -RL ${candidates} ${prefix}/candidates

    varnet-predict \\
        --sample_name ${prefix} \\
        --normal_bam ${input_normal} \\
        --tumor_bam ${input_tumor} \\
        --reference ${fasta} \\
        --output_dir . \\
        --processes ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/predictions/snvs ${prefix}/predictions/indels
    touch ${prefix}/predictions/snvs/Predictions.csv
    touch ${prefix}/predictions/indels/Predictions.csv
    echo "" | gzip > ${prefix}/${prefix}.vcf.gz
    """
}
