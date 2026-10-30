// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process VARNET_PREDICT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/varnet:1.5.4--pyhdfd78af_0':
        'quay.io/biocontainers/varnet:1.5.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(input_tumor), path(index_tumor), path(input_normal), path(index_normal), path(candidates, stageAs: 'candidates')
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    // tumor-only panel of normals: varnet-predict finds it by this name in --resource_dir
    tuple val(meta4), path(pon, stageAs: 'varnet_resources/1000g_pon.hg38.vcf.gz')

    output:
    tuple val(meta), path("${prefix}/${prefix}.vcf.gz"), emit: vcf
    tuple val(meta), path("${prefix}/predictions"), emit: predictions
    tuple val("${task.process}"), val("varnet"), eval("varnet-predict --version | sed 's/VarNet //'"), emit: versions_varnet, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def normal = input_normal ? "--normal_bam ${input_normal}" : ""
    def resources = pon ? "--resource_dir varnet_resources" : ""
    if (!input_normal && !resources) {
        error "VARNET_PREDICT tumor-only mode needs the panel of normals from https://doi.org/10.6084/m9.figshare.34420128 (or varnet-download-resources)."
    }
    """
    export TF_CPP_MIN_LOG_LEVEL=3

    # varnet-predict expects the candidates from varnet-filter at <output_dir>/<sample_name>/candidates
    mkdir -p ${prefix}
    cp -RL ${candidates} ${prefix}/candidates

    varnet-predict \\
        --sample_name ${prefix} \\
        ${normal} \\
        --tumor_bam ${input_tumor} \\
        --reference ${fasta} \\
        --output_dir . \\
        --processes ${task.cpus} \\
        ${resources} \\
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
