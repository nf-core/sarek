process VARLOCIRAPTOR_CALLVARIANTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9a/9ac0825c21b2cbaf9535ffe443e53a0bb4d61596cafcb5a5b444dfb31b945ab2/data'
        : 'community.wave.seqera.io/library/varlociraptor:8.9.3--fa2ce5da2782669c'}"

    input:
    tuple val(meta), path(vcfs), path(scenario), val(scenario_aliases)

    output:
    tuple val(meta), path("*.bcf"), emit: bcf
    tuple val("${task.process}"), val('varlociraptor'), eval("varlociraptor --version | sed 's/^varlociraptor //'"), topic: versions, emit: versions_varlociraptor

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_called"

    //If we use a scenario file and if there is more than 1 vcf, then collect scenario_aliaes and vcfs to scenario_alias_0=vcf_0 scenario_alias_1=vcf_1, etc
    //If we use a scenario file and if there is exactly 1 vcf, then scenario_alias=vcf
    def scenario_samples = vcfs instanceof List && vcfs.size() > 1 ? [scenario_aliases, vcfs].transpose().collect { files -> "${files[0]}=${files[1]}" }.join(' ') : "${scenario_aliases}=${vcfs}"
    """
    varlociraptor call variants \\
        --output ${prefix}.bcf \\
        generic --scenario ${scenario} --obs ${scenario_samples} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_called"
    """
    touch ${prefix}.bcf
    """
}
