process SNPSIFT_ANNMEM {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/35/3570381a93c22449d48bdaa85097c5e8a075e90437565546acb2e40a29171bca/data'
        : 'community.wave.seqera.io/library/snpsift:5.3.0a--67d3871d6f67ac2b'}"

    input:
    tuple val(meta), path(vcf), path(vcf_tbi)
    path databases                                  // List of annotation VCF files
    path dbs_tbis                                   // List of corresponding index files
    path db_vardbs, stageAs: 'db_vardbs/*'         // List of .snpsift.vardb directories
    val db_configs                                  // List of maps: [[fields: 'field1,field2', prefix: 'DB1_'], ...]

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def db_list = databases instanceof List ? databases : [databases]

    // Build -dbfile arguments for each database
    def dbfile_args = []
    db_list.eachWithIndex { db, idx ->
        def config = db_configs && db_configs[idx] ? db_configs[idx] : [:]
        def db_arg = "-dbfile ${db}"
        if (config.fields) {
            db_arg += " -fields ${config.fields}"
        }
        if (config.prefix) {
            db_arg += " -prefix ${config.prefix}"
        }
        dbfile_args << db_arg
    }

    """
    SnpSift \\
        annmem \\
        ${args} \\
        ${dbfile_args.join(' \\\n        ')} \\
        ${vcf} \\
        | bgzip -c > ${prefix}.vcf.gz

    tabix -p vcf ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpsift: \$( echo \$(SnpSift split -h 2>&1) | sed 's/^.*version //' | sed 's/(.*//' | sed 's/t//g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpsift: \$( echo \$(SnpSift split -h 2>&1) | sed 's/^.*version //' | sed 's/(.*//' | sed 's/t//g' )
    END_VERSIONS
    """
}
