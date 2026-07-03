process ENSEMBLVEP_DOWNLOAD {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/edd02dfaf968d06c808e3c208d5b3e86afb4259590bfa6e5499965ef3bc81881/data'
        : 'community.wave.seqera.io/library/ensembl-vep_perl-math-cdf_htslib:efd9a6d1c5f218a9'}"

    input:
    tuple val(meta), val(assembly), val(species), val(cache_version)
    val(preflight_check)

    output:
    tuple val(meta), path(prefix), emit: cache
    tuple val("${task.process}"), val('ensemblvep'), eval("vep --help | sed -n '/ensembl-vep/s/.*: //p'"), topic: versions, emit: versions_ensemblvep
    tuple val("${task.process}"), val('perl-math-cdf'), eval("perl -MMath::CDF -e 'print \$Math::CDF::VERSION'"), topic: versions, emit: versions_perlmathcdf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: 'vep_cache'
    def filename = "${species}_vep_${cache_version}_${assembly}.tar.gz"
    def checksums_url = "https://ftp.ensembl.org/pub/release-${cache_version}/variation/indexed_vep_cache/CHECKSUMS"
    """
    if [ "${preflight_check}" = "true" ]; then
        perl -MHTTP::Tiny -e '
            my \$r = HTTP::Tiny->new(timeout => 30)->get("${checksums_url}");
            \$r->{success} or die "Failed to fetch CHECKSUMS (HTTP \$r->{status})\\n";
            \$r->{content} =~ /\\Q${filename}\\E/ or die "${filename} not found in CHECKSUMS\\n";
            print "Pre-flight OK: ${filename} found in CHECKSUMS\\n";
        '
    fi

    vep_install \\
        --CACHEDIR ${prefix} \\
        --SPECIES ${species} \\
        --ASSEMBLY ${assembly} \\
        --CACHE_VERSION ${cache_version} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: 'vep_cache'
    """
    mkdir ${prefix}
    """
}
