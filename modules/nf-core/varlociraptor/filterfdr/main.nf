process VARLOCIRAPTOR_FILTERFDR {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/4788613cdc4faf20c5ebb24da24575ec3d80a8a52a0daf1c986b742444fe6c99/data'
:         'community.wave.seqera.io/library/varlociraptor:8.12.0--294eead5eab5a1dc' }"

    input:
    tuple val(meta), path(vcf), val(events), val(fdr)

    output:
    tuple val(meta), path("*.bcf"), emit: bcf
    tuple val("${task.process}"), val('varlociraptor'), eval("varlociraptor --version | sed 's/^varlociraptor //'"), topic: versions, emit: versions_varlociraptor

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.fdr-controlled"
    def mode = args.contains("--mode global-smart") ? "global-smart" : "local-smart"
    def args_corrected = args.replace('--mode global-smart', '').replace('--mode local-smart', '').trim()
    """
    varlociraptor filter-calls control-fdr \\
        --mode ${mode} \\
        ${vcf} \\
        --events ${events} \\
        --fdr ${fdr} \\
        ${args_corrected} \\
        > ${prefix}.bcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.fdr-controlled"
    """
    touch ${prefix}.bcf
    """
}
