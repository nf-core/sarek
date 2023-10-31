//
// INITIALIZE ANNOTATION CACHE
//

// Initialize channels based on params or indices that were just built
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
// Condition is based on params.step and params.tools
// If and extra condition exists, it's specified in comments

workflow INITIALIZE_ANNOTATION_CACHE {
    take:
    snpeff_enabled
    snpeff_cache
    snpeff_genome
    snpeff_db
    vep_enabled
    vep_cache
    vep_species
    vep_cache_version
    vep_genome
    help_message

    main:
    if (snpeff_enabled) {
        def snpeff_annotation_cache_key = (snpeff_cache == "s3://annotation-cache/snpeff_cache/") ? "${snpeff_genome}.${snpeff_db}/" : ""
        def snpeff_cache_dir =  "${snpeff_annotation_cache_key}${snpeff_genome}.${snpeff_db}"
        def snpeff_cache_path_full = file("$snpeff_cache/$snpeff_cache_dir", type: 'dir')
        if ( !snpeff_cache_path_full.exists() || !snpeff_cache_path_full.isDirectory() ) {
            if (snpeff_cache == "s3://annotation-cache/snpeff_cache/") {
                error("This path is not available within annotation-cache.\nPlease check https://annotation-cache.github.io/ to create a request for it.")
            } else {
                error("Path provided with SnpEff cache is invalid.\nMake sure there is a directory named ${snpeff_cache_dir} in ${snpeff_cache}./n${help_message}")
            }
        }
        snpeff_cache = Channel.fromPath(file("${snpeff_cache}/${snpeff_annotation_cache_key}"), checkIfExists: true).collect()
            .map{ cache -> [ [ id:"${snpeff_genome}.${snpeff_db}" ], cache ] }
    } else snpeff_cache = []

    if (vep_enabled) {
        def vep_annotation_cache_key = (vep_cache == "s3://annotation-cache/vep_cache/") ? "${vep_cache_version}_${vep_genome}/" : ""
        def vep_cache_dir = "${vep_annotation_cache_key}${vep_species}/${vep_cache_version}_${vep_genome}"
        def vep_cache_path_full = file("$vep_cache/$vep_cache_dir", type: 'dir')
        if ( !vep_cache_path_full.exists() || !vep_cache_path_full.isDirectory() ) {
            if (vep_cache == "s3://annotation-cache/vep_cache/") {
                error("This path is not available within annotation-cache.\nPlease check https://annotation-cache.github.io/ to create a request for it.")
            } else {
                error("Path provided with VEP cache is invalid.\nMake sure there is a directory named ${vep_cache_dir} in ${vep_cache}./n${help_message}")
            }
        }
        ensemblvep_cache = Channel.fromPath(file("${vep_cache}/${vep_annotation_cache_key}"), checkIfExists: true).collect()
    } else ensemblvep_cache = []

    emit:
    ensemblvep_cache // channel: [ meta, cache ]
    snpeff_cache     // channel: [ meta, cache ]
}
