process BWA_INDEX {
    tag "$fasta"
    // NOTE requires 5.37N memory where N is the size of the database
    // source: https://bio-bwa.sourceforge.net/bwa.shtml#8
    memory { 7.B * fasta.size() }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d7/d7e24dc1e4d93ca4d3a76a78d4c834a7be3985b0e1e56fddd61662e047863a8a/data' :
        'community.wave.seqera.io/library/bwa_htslib_samtools:83b50ff84ead50d0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bwa"), emit: index
    tuple val("${task.process}"), val('bwa'), eval('bwa 2>&1 | sed -n "s/^Version: //p"'), topic: versions, emit: versions_bwa

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    def args   = task.ext.args ?: ''
    """
    mkdir bwa
    bwa \\
        index \\
        $args \\
        -p bwa/${prefix} \\
        $fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    mkdir bwa
    touch bwa/${prefix}.amb
    touch bwa/${prefix}.ann
    touch bwa/${prefix}.bwt
    touch bwa/${prefix}.pac
    touch bwa/${prefix}.sa
    """
}
