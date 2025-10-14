process ENSEMBLVEP_VEP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4b/4b5a8c173dc9beaa93effec76b99687fc926b1bd7be47df5d6ce19d7d6b4d6b7/data'
        : 'community.wave.seqera.io/library/ensembl-vep:115.2--90ec797ecb088e9a'}"

    input:
    tuple val(meta), path(vcf), path(custom_extra_files)
    val genome
    val species
    val cache_version
    path cache
    tuple val(meta2), path(fasta)
    path extra_files

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf, optional: true
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi, optional: true
    tuple val(meta), path("*.tab.gz"), emit: tab, optional: true
    tuple val(meta), path("*.json.gz"), emit: json, optional: true
    path "*.html", emit: report, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def file_extension = args.contains("--vcf") ? 'vcf' : args.contains("--json") ? 'json' : args.contains("--tab") ? 'tab' : 'vcf'
    def compress_cmd = args.contains("--compress_output") ? '' : '--compress_output bgzip'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dir_cache = cache ? "\${PWD}/${cache}" : "/.vep"
    def reference = fasta ? "--fasta ${fasta}" : ""
    def create_index = file_extension == "vcf" ? "tabix ${args2} ${prefix}.${file_extension}.gz" : ""
    """
    vep \\
        -i ${vcf} \\
        -o ${prefix}.${file_extension}.gz \\
        ${args} \\
        ${compress_cmd} \\
        ${reference} \\
        --assembly ${genome} \\
        --species ${species} \\
        --cache \\
        --cache_version ${cache_version} \\
        --dir_cache ${dir_cache} \\
        --fork ${task.cpus}

    ${create_index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_extension = args.contains("--vcf") ? 'vcf' : args.contains("--json") ? 'json' : args.contains("--tab") ? 'tab' : 'vcf'
    def create_index = file_extension == "vcf" ? "touch ${prefix}.${file_extension}.gz.tbi" : ""
    """
    echo "" | gzip > ${prefix}.${file_extension}.gz
    ${create_index}
    touch ${prefix}_summary.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
