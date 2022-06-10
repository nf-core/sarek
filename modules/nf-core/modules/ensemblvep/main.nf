process ENSEMBLVEP {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::ensembl-vep=104.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ensembl-vep:104.3--pl5262h4a94de4_0' :
        'quay.io/biocontainers/ensembl-vep:104.3--pl5262h4a94de4_0' }"

    input:
    tuple val(meta), path(vcf)
    val   vep_output
    val   genome
    val   species
    val   cache_version
    path  cache
    path  extra_files

    output:
    tuple val(meta), path("*.ann.vcf")  , optional:true, emit: vcf
    tuple val(meta), path("*.ann.tab")  , optional:true, emit: tab
    tuple val(meta), path("*.ann.json") , optional:true, emit: json
    path "*.summary.html"               , emit: report
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def suffix = task.ext.suffix ?: ''
    def prefix = args.prefix ? "$args.prefix" : "${meta.id}"
    def dir_cache = cache ? "\${PWD}/${cache}" : "/.vep"
    """
    mkdir $prefix

    vep \\
        -i $vcf \\
        -o ${prefix}${suffix}.ann.$vep_output \\
        $args \\
        --assembly $genome \\
        --species $species \\
        --cache \\
        --cache_version $cache_version \\
        --dir_cache $dir_cache \\
        --fork $task.cpus \\
        --stats_file ${prefix}.summary.html \\
        --$vep_output

    rm -rf $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}
