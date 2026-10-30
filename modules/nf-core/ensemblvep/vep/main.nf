process ENSEMBLVEP_VEP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6f/6febf703e6df28e6259f8b382a75e16cd7039be1da722a8ee716004479b82458/data'
        : 'community.wave.seqera.io/library/ensembl-vep_perl-math-cdf_htslib_gzip_tar:43f77a48cc79bf38'}"

    input:
    tuple val(meta), path(vcf), path(custom_extra_files)
    val genome
    val species
    val cache_version
    tuple val(meta2), path(cache)
    tuple val(meta3), path(fasta)
    path extra_files
    tuple path(gtf), path(gtf_tbi) // optional: [ path(gtf), path(gtf_tbi) ] -- mutually exclusive with cache

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf, optional: true
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: tbi, optional: true
    tuple val(meta), path("${prefix}.tab.gz"), emit: tab, optional: true
    tuple val(meta), path("${prefix}.json.gz"), emit: json, optional: true
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
    prefix = task.ext.prefix ?: "${meta.id}"
    def dir_cache = cache ? "\${PWD}/${cache}" : "/.vep"
    def reference = fasta ? "--fasta ${fasta}" : ""
    // gtf and cache are mutually exclusive annotation sources: use a GTF
    // (e.g. for a custom/non-standard reference with no Ensembl cache)
    // when supplied, otherwise fall back to the existing cache behaviour
    // unchanged.
    def annotation_source = gtf ? "--gtf ${gtf}" : "--cache --cache_version ${cache_version} --dir_cache ${dir_cache}"
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
        ${annotation_source} \\
        --fork ${task.cpus}

    ${create_index}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def file_extension = args.contains("--vcf") ? 'vcf' : args.contains("--json") ? 'json' : args.contains("--tab") ? 'tab' : 'vcf'
    def create_index = file_extension == "vcf" ? "touch ${prefix}.${file_extension}.gz.tbi" : ""
    """
    echo "" | gzip > ${prefix}.${file_extension}.gz
    ${create_index}
    touch ${prefix}_summary.html
    """
}
