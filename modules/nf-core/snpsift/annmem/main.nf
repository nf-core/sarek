process SNPSIFT_ANNMEM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a1/a116bb44e388ca83fea78d82fe8bdfd5cf3557254e2ec7dd3f1f17354880638c/data' :
        'community.wave.seqera.io/library/htslib_snpsift:ace461dff1cfc121' }"

    input:
    tuple val(meta), path(vcf), path(vcf_tbi)
    tuple path(db_vcf), path(db_vcf_tbi), path(db_vardb), val(db_fields), val(db_prefixes)

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('snpsift'), eval("SnpSift -version 2>&1 | grep -oE '[0-9]+\\.[0-9]+[a-z]?'"), emit: versions_snpsift, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def dbs      = db_vcf instanceof List ? db_vcf : [db_vcf]
    def all_fields = db_fields instanceof List ? db_fields : [db_fields]
    def prefixes = db_prefixes instanceof List ? db_prefixes : [db_prefixes]

    // db_vardb is staged as input so it's present in the work directory;
    // SnpSift finds it automatically next to the VCF as {vcf}.snpsift.vardb/
    def dbfile_args = dbs.withIndex().collect { db, i ->
        def f = all_fields[i]
        def fields = f instanceof List ? f.join(',') : f
        def p = prefixes[i]
        "-dbfile ${db}${fields ? " -fields ${fields}" : ''}${p ? " -prefix ${p}" : ''}"
    }

    """
    SnpSift \\
        annmem \\
        ${args} \\
        ${dbfile_args.join(' \\\n        ')} \\
        ${vcf} \\
        | bgzip -c > ${prefix}.vcf.gz

    tabix -p vcf ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | bgzip -c > ${prefix}.vcf.gz
    echo "" | gzip > ${prefix}.vcf.gz.tbi
    """
}
