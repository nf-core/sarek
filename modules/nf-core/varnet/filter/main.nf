// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process VARNET_FILTER {
    tag "$meta.id"
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/varnet:1.5.4--pyhdfd78af_0':
        'quay.io/biocontainers/varnet:1.5.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(input_tumor), path(index_tumor), path(input_normal), path(index_normal)
    tuple val(meta2), path(intervals)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)
    // tumor-only germline filtering resources: varnet-filter finds them by these names in --resource_dir
    tuple val(meta5), path(dbsnp, stageAs: 'varnet_resources/GrCh38-dbSNP-build-156-common-variants-0-indexed-CHROMPOS.npz')
    tuple val(meta6), path(gnomad, stageAs: 'varnet_resources/combined_gnomAD_v4.1_genomes_exomes_dict.hg38.chunks.npz')

    output:
    tuple val(meta), path("${prefix}/candidates"), emit: candidates
    tuple val("${task.process}"), val("varnet"), eval("varnet-filter --version | sed 's/VarNet //'"), emit: versions_varnet, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def regions = intervals ? "--region_bed ${intervals}" : ""
    def normal = input_normal ? "--normal_bam ${input_normal}" : ""
    def resources = dbsnp && gnomad ? "--resource_dir varnet_resources" : ""
    if (!input_normal && !resources) {
        error "VARNET_FILTER tumor-only mode needs the dbSNP and gnomAD germline filtering resources from https://doi.org/10.6084/m9.figshare.34420128 (or varnet-download-resources)."
    }
    """
    export TF_CPP_MIN_LOG_LEVEL=3

    varnet-filter \\
        --sample_name ${prefix} \\
        ${normal} \\
        --tumor_bam ${input_tumor} \\
        --reference ${fasta} \\
        --output_dir . \\
        --processes ${task.cpus} \\
        ${regions} \\
        ${resources} \\
        ${args}

    # varnet-filter joins its batch files in directory listing order, which differs between
    # filesystems, so sort the candidates to make the output reproducible
    # (-snv or -indel in ext.args writes only one of the two files)
    for candidates in ${prefix}/candidates/*/Positions.csv; do
        LC_ALL=C sort -k1,1 -k2,2n -o \$candidates \$candidates
    done
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/candidates/snvs ${prefix}/candidates/indels
    touch ${prefix}/candidates/snvs/Positions.csv
    touch ${prefix}/candidates/indels/Positions.csv
    """
}
