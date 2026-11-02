//
// DOWNLOAD CACHE SNPEFF VEP
//

include { ENSEMBLVEP_DOWNLOAD } from '../../../modules/nf-core/ensemblvep/download'
include { SNPEFF_DOWNLOAD     } from '../../../modules/nf-core/snpeff/download'

workflow CACHE_DOWNLOAD_ENSEMBLVEP_SNPEFF {
    take:
    ensemblvep_info
    snpeff_info
    ensemblvep_preflight_check

    main:
    ENSEMBLVEP_DOWNLOAD(ensemblvep_info, ensemblvep_preflight_check)
    SNPEFF_DOWNLOAD(snpeff_info)

    emit:
    ensemblvep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.collect() // channel: [ meta, cache ]
    snpeff_cache     = SNPEFF_DOWNLOAD.out.cache.collect() // channel: [ meta, cache ]
}
