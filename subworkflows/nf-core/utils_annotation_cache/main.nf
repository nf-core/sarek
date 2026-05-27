//
// Subworkflow to help using annotation cache
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow UTILS_ANNOTATION_CACHE {
    take:
    ensemblvep_cache
    ensemblvep_cache_version
    ensemblvep_custom_args
    ensemblvep_genome
    ensemblvep_species
    ensemblvep_enabled
    snpeff_cache
    snpeff_db
    snpeff_enabled
    help_message

    main:
    if (ensemblvep_enabled) {
        def ensemblvep_annotation_cache_key = isCloudUrl(ensemblvep_cache) ? "${ensemblvep_cache_version}_${ensemblvep_genome}/" : ""
        def ensemblvep_species_suffix = ensemblvep_custom_args.contains("--merged") ? '_merged' : (ensemblvep_custom_args.contains("--refseq") ? '_refseq' : '')
        def ensemblvep_cache_dir = "${ensemblvep_annotation_cache_key}${ensemblvep_species}${ensemblvep_species_suffix}/${ensemblvep_cache_version}_${ensemblvep_genome}"
        def ensemblvep_cache_path_full = file("${ensemblvep_cache}/${ensemblvep_cache_dir}", type: 'dir')
        if (!ensemblvep_cache_path_full.exists() || !ensemblvep_cache_path_full.isDirectory()) {
            if (ensemblvep_cache == "s3://annotation-cache/vep_cache/") {
                error("This path is not available within annotation-cache.\nPlease check https://annotation-cache.github.io/ to create a request for it.")
            }
            else {
                error("Path provided with ENSEMBLVEP cache is invalid.\nMake sure there is a directory named ${ensemblvep_cache_dir} in ${ensemblvep_cache}.\n${help_message}")
            }
        }
        ensemblvep_cache = channel.fromPath(file("${ensemblvep_cache}/${ensemblvep_annotation_cache_key}"), checkIfExists: true)
            .collect()
            .map { cache -> [[id: "${ensemblvep_cache_version}_${ensemblvep_genome}"], cache] }
    }
    else {
        ensemblvep_cache = []
    }

    if (snpeff_enabled) {
        def snpeff_annotation_cache_key = isCloudUrl(snpeff_cache) ? "${snpeff_db}/" : ""
        def snpeff_cache_dir = "${snpeff_annotation_cache_key}${snpeff_db}"
        def snpeff_cache_path_full = file("${snpeff_cache}/${snpeff_cache_dir}", type: 'dir')
        if (!snpeff_cache_path_full.exists() || !snpeff_cache_path_full.isDirectory()) {
            if (snpeff_cache == "s3://annotation-cache/snpeff_cache/") {
                error("This path is not available within annotation-cache.\nPlease check https://annotation-cache.github.io/ to create a request for it.")
            }
            else {
                error("Path provided with SnpEff cache is invalid.\nMake sure there is a directory named ${snpeff_cache_dir} in ${snpeff_cache}.\n${help_message}")
            }
        }
        snpeff_cache = channel.fromPath(file("${snpeff_cache}/${snpeff_annotation_cache_key}"), checkIfExists: true)
            .collect()
            .map { cache -> [[id: "${snpeff_db}"], cache] }
    }
    else {
        snpeff_cache = []
    }

    emit:
    ensemblvep_cache // channel: [ meta, cache ]
    snpeff_cache // channel: [ meta, cache ]
}

// Helper function to check if cache path is from any cloud provider
def isCloudUrl(cache_url) {
    return cache_url.startsWith("s3://") || cache_url.startsWith("gs://") || cache_url.startsWith("az://")
}
