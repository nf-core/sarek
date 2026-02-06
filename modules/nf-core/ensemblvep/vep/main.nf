process ENSEMBLVEP_VEP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3d/3da6e21cbf9803529421d7e136d1ebec5ff71ec50e0d996eda2ce11ec2c19bf9/data'
        : 'community.wave.seqera.io/library/ensembl-vep_perl-math-cdf:1e13f65f931a6954'}"

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
    tuple val(meta), val("${task.process}"), val('ensemblvep'), path("*.html"), topic: multiqc_files, emit: report, optional: true
    tuple val("${task.process}"), val('ensemblvep'), eval("vep --help | sed -n '/ensembl-vep/s/.*: //p'"), topic: versions, emit: versions_ensemblvep
    tuple val("${task.process}"), val('tabix'), eval("tabix -h 2>&1 | grep -oP 'Version:\\s*\\K[^\\s]+'"), topic: versions, emit: versions_tabix
    tuple val("${task.process}"), val('perl-math-cdf'), eval("perl -MMath::CDF -e 'print \\\$Math::CDF::VERSION'"), topic: versions, emit: versions_perlmathcdf

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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_extension = args.contains("--vcf") ? 'vcf' : args.contains("--json") ? 'json' : args.contains("--tab") ? 'tab' : 'vcf'
    def create_index = file_extension == "vcf" ? "touch ${prefix}.${file_extension}.gz.tbi" : ""
    """
    echo "" | gzip > ${prefix}.${file_extension}.gz
    ${create_index}
    touch ${prefix}_summary.html
    """
}
