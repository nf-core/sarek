process SNPSIFT_ANNMEM {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/35/3570381a93c22449d48bdaa85097c5e8a075e90437565546acb2e40a29171bca/data'
        : 'community.wave.seqera.io/library/htslib_snpsift:4df051c7ff79f7f9'}"

    input:
    tuple val(meta), path(vcf), path(vcf_tbi)
    tuple path(databases), path(dbs_tbis), path(db_vardbs, stageAs: 'db_vardbs/*'), val(db_fields), val(db_prefixes)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    tuple val("${task.process}"), val('snpsift'), eval("SnpSift split -h 2>&1 | head -1 | sed 's/^.*version //' | sed 's/(.*//' | sed 's/t//g'"), topic: 'versions'

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def db_list = databases instanceof List ? databases : [databases]
    def fields_list = db_fields instanceof List ? db_fields : [db_fields]
    def prefix_list = db_prefixes instanceof List ? db_prefixes : [db_prefixes]

    // Build -dbfile arguments for each database (semicolons converted to commas for SnpSift)
    def dbfile_args = db_list.withIndex().collect { db, idx ->
        def fields = fields_list[idx]?.replace(';', ',')
        def db_prefix = prefix_list[idx]
        "-dbfile ${db}" +
            (fields ? " -fields ${fields}" : '') +
            (db_prefix ? " -prefix ${db_prefix}" : '')
    }

    """
    # Link database directories from db_vardbs/ to current directory
    # so SnpSift can find them next to the VCF files
    if [ -d "db_vardbs" ]; then
        for vardb in db_vardbs/*.snpsift.vardb; do
            if [ -d "\$vardb" ]; then
                ln -s "\$vardb" .
            fi
        done
    fi

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
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
