//
// ANNOTATION CACHE INITIALISATION
//

// Initialise channels based on params or indices that were just built
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
// Condition is based on params.step and params.tools
// If and extra condition exists, it's specified in comments

workflow ANNOTATION_CACHE_INITIALISATION {
    take:
    snpeff_enabled
    snpeff_cache_in
    snpeff_db
    vep_enabled
    vep_cache
    vep_species
    vep_cache_version
    vep_genome
    vep_custom_args
    help_message

    main:
    if (snpeff_enabled) {
        snpeff_cache = snpeff_db.map { _id, snpeff_db_ ->
            def snpeff_annotation_cache_key = snpeff_cache_in == "s3://annotation-cache/snpeff_cache/" ? "${snpeff_db_}/" : ""
            def snpeff_cache_dir = "${snpeff_annotation_cache_key}${snpeff_db_}"
            def snpeff_cache_path_full = file("${snpeff_cache_in}/${snpeff_cache_dir}", type: 'dir')
            if (!snpeff_cache_path_full.exists() || !snpeff_cache_path_full.isDirectory()) {
                if (snpeff_cache_in == "s3://annotation-cache/snpeff_cache/") {
                    error("This path is not available within annotation-cache.\nPlease check https://annotation-cache.github.io/ to create a request for it.")
                }
                else {
                    error("Path provided with SnpEff cache is invalid.\nMake sure there is a directory named ${snpeff_cache_dir} in ${snpeff_cache_in}./n${help_message}")
                }
            }
            [[id: snpeff_db_], file("${snpeff_cache_in}/${snpeff_annotation_cache_key}")]
        }
    }
    else {
        snpeff_cache = []
    }

    if (vep_enabled) {
        ensemblvep_cache = vep_cache_version
            .join(vep_species)
            .join(vep_genome)
            .groupTuple()
            .map { _id, vep_cache_version_, vep_species_, vep_genome_ ->
                def vep_annotation_cache_key = vep_cache == "s3://annotation-cache/vep_cache/" ? "${vep_cache_version_[0]}_${vep_genome_[0]}/" : ""
                def vep_species_suffix = vep_custom_args.contains("--merged") ? '_merged' : (vep_custom_args.contains("--refseq") ? '_refseq' : '')
                def vep_cache_dir = "${vep_annotation_cache_key}${vep_species_[0]}${vep_species_suffix}/${vep_cache_version_[0]}_${vep_genome_[0]}"
                def vep_cache_path_full = file("${vep_cache}/${vep_cache_dir}", type: 'dir')

                if (!vep_cache_path_full.exists() || !vep_cache_path_full.isDirectory()) {
                    if (vep_cache == "s3://annotation-cache/vep_cache/") {
                        error("This path is not available within annotation-cache.\nPlease check https://annotation-cache.github.io/ to create a request for it.")
                    }
                    else {
                        error("Path provided with VEP cache is invalid.\nMake sure there is a directory named ${vep_cache_dir} in ${vep_cache}./n${help_message}")
                    }
                }
                [file("${vep_cache}/${vep_annotation_cache_key}")]
            }
    }
    else {
        ensemblvep_cache = []
    }

    emit:
    ensemblvep_cache // channel: [ meta, cache ]
    snpeff_cache     // channel: [ meta, cache ]
}
