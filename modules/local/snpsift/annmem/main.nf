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
    val(create)

    output:
    tuple val(meta), path("*.snpsift.vardb"), emit: database, optional: true
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf, optional: true
    tuple val("${task.process}"), val('snpsift'), eval("SnpSift split -h 2>&1 | head -1 | sed 's/^.*version //' | sed 's/(.*//' | sed 's/t//g'"), topic: 'versions'

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (create) {
        def db = db_vcf instanceof List ? db_vcf[0] : db_vcf
        def f = db_fields instanceof List ? db_fields[0] : db_fields
        def fields = f instanceof List ? f.join(',') : f?.replace(';', ',')
        """
        SnpSift \\
            annmem \\
            -create \\
            ${args} \\
            -dbfile ${db} \\
            ${fields ? "-fields ${fields}" : ""}
        """
    } else {
        def dbs = db_vcf instanceof List ? db_vcf : [db_vcf]
        def vardbs = db_vardb instanceof List ? db_vardb : [db_vardb]
        def all_fields = db_fields instanceof List ? db_fields : [db_fields]
        def prefixes = db_prefixes instanceof List ? db_prefixes : [db_prefixes]

        def dbfile_args = dbs.withIndex().collect { db, i ->
            def dbfile = vardbs[i] ?: db
            def f = all_fields[i]
            def fields = f instanceof List ? f.join(',') : f?.replace(';', ',')
            def p = prefixes[i]
            "-dbfile ${dbfile}${fields ? " -fields ${fields}" : ''}${p ? " -prefix ${p}" : ''}"
        }

        """
        SnpSift \\
            annmem \\
            ${args} \\
            ${dbfile_args.join(' \\\n            ')} \\
            ${vcf} \\
            | bgzip -c > ${prefix}.vcf.gz

        tabix -p vcf ${prefix}.vcf.gz
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def db = db_vcf instanceof List ? db_vcf[0] : db_vcf
    if (create) {
        """
        mkdir -p ${db}.snpsift.vardb
        touch ${db}.snpsift.vardb/chr1.snpsift.df
        """
    } else {
        """
        touch ${prefix}.vcf.gz
        touch ${prefix}.vcf.gz.tbi
        """
    }
}
