process SNPSIFT_ANNMEM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a1/a116bb44e388ca83fea78d82fe8bdfd5cf3557254e2ec7dd3f1f17354880638c/data' :
        'community.wave.seqera.io/library/htslib_snpsift:ace461dff1cfc121' }"

    input:
    tuple val(meta), path(vcf), path(vcf_tbi), val(fields)
    tuple path(databases), path(dbs_tbis), path(db_vardbs, stageAs: 'db_vardbs/*'), val(db_fields), val(db_prefixes)

    output:
    tuple val(meta), path("*.snpsift.vardb"), emit: database, optional: true
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf, optional: true
    tuple val("${task.process}"), val('snpsift'), eval("SnpSift split -h 2>&1 | head -1 | sed 's/^.*version //' | sed 's/(.*//' | sed 's/t//g'"), topic: 'versions'

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def create_mode = task.ext.create ?: false

    if (create_mode) {
        // Database creation mode
        // Convert semicolons to commas for SnpSift (semicolons used in CSV to avoid delimiter conflict)
        def fields_csv = fields ? fields.toString().replace(';', ',') : ''
        def fields_arg = fields_csv ? "-fields ${fields_csv}" : ""
        """
        SnpSift \\
            annmem \\
            -create \\
            ${args} \\
            -dbfile ${vcf} \\
            ${fields_arg}
        """
    } else {
        // Annotation mode
        def db_list = databases instanceof List ? databases : [databases]
        def fields_list = db_fields instanceof List ? db_fields : [db_fields]
        def prefix_list = db_prefixes instanceof List ? db_prefixes : [db_prefixes]

        // Build -dbfile arguments for each database (semicolons converted to commas for SnpSift)
        def dbfile_args = db_list.withIndex().collect { db, idx ->
            def db_fields_str = fields_list[idx]?.toString()?.replace(';', ',')
            def db_prefix = prefix_list[idx]
            "-dbfile ${db}" +
                (db_fields_str ? " -fields ${db_fields_str}" : '') +
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
            ${dbfile_args.join(' \\\n            ')} \\
            ${vcf} \\
            | bgzip -c > ${prefix}.vcf.gz

        tabix -p vcf ${prefix}.vcf.gz
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def create_mode = task.ext.create ?: false
    if (create_mode) {
        """
        mkdir -p ${vcf}.snpsift.vardb
        touch ${vcf}.snpsift.vardb/chr1.snpsift.df
        """
    } else {
        """
        touch ${prefix}.vcf.gz
        touch ${prefix}.vcf.gz.tbi
        """
    }
}
