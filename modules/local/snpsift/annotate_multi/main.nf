process SNPSIFT_ANNOTATE_MULTI {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/35/3570381a93c22449d48bdaa85097c5e8a075e90437565546acb2e40a29171bca/data'
        : 'community.wave.seqera.io/library/snpsift:5.3.0a--67d3871d6f67ac2b'}"

    input:
    tuple val(meta), path(vcf), path(vcf_tbi)
    val db_list                                     // List of maps: [[name:'dbsnp', file:'...', fields:'', id:true], ...]
    path databases                                  // All database VCF files
    path dbs_tbis                                   // All database TBI files

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Build pipeline of SnpSift commands
    def snpsift_cmds = []
    db_list.each { db ->
        def db_args = []
        if (db.id) db_args << '-id'
        if (db.fields) db_args << "-info ${db.fields}"

        snpsift_cmds << "SnpSift annotate ${db_args.join(' ')} ${db.file}"
    }

    """
    # Chain all SnpSift commands via pipes
    zcat ${vcf} | \\
        ${snpsift_cmds.join(' | \\\n        ')} | \\
        bgzip -c > ${prefix}.vcf.gz

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
