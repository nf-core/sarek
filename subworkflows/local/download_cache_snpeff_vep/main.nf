//
// DOWNLOAD CACHE SNPEFF VEP
//

// Initialize channels based on params or indices that were just built
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
// Condition is based on params.step and params.tools
// If and extra condition exists, it's specified in comments

include { ENSEMBLVEP_DOWNLOAD } from '../../../modules/nf-core/ensemblvep/download/main'
include { SNPEFF_DOWNLOAD     } from '../../../modules/nf-core/snpeff/download/main'

workflow DOWNLOAD_CACHE_SNPEFF_VEP {
    take:
    ensemblvep_info
    snpeff_info

    main:
    versions = Channel.empty()

    ENSEMBLVEP_DOWNLOAD(ensemblvep_info)
    SNPEFF_DOWNLOAD(snpeff_info)

    // Gather versions of all tools used
    versions = versions.mix(ENSEMBLVEP_DOWNLOAD.out.versions)
    versions = versions.mix(SNPEFF_DOWNLOAD.out.versions)

    emit:
    ensemblvep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.collect()  // channel: [ meta, cache ]
    snpeff_cache     = SNPEFF_DOWNLOAD.out.cache.collect()      // channel: [ meta, cache ]

    versions // channel: [ versions.yml ]
}
