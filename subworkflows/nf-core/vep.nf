/*
 * Run VEP to annotate VCF files
 */

params.vep_options     = [:]
params.bgziptabix_vep  = [:]

include { VEP }              from '../../modules/nf-core/software/vep/main'              addParams(options: params.vep_options)
include { TABIX_BGZIPTABIX } from '../../modules/nf-core/software/tabix/bgziptabix/main' addParams(options: params.bgziptabix_vep)

workflow VEP_ANNOTATE {
    take:
    vcf               // channel: [ val(meta), vcf, tbi ]
    vep_genome        //   value: which genome
    vep_species       //   value: which species
    vep_cache_version //   value: which cache version
    use_cache         //    bool: use cache (default: false)
    vep_cache         //    path: path_to_vep_cache (optionnal)
    vep_tag           //   value: tag of docker image (optionnal)
    // skip   // boolean: true/false

    main:
    params.use_cache = false
    use_cache = params.use_cache
    params.vep_tag = false
    vep_tag = params.vep_tag

    VEP(vcf, vep_genome, vep_species, vep_cache_version, use_cache, vep_cache, vep_tag)
    TABIX_BGZIPTABIX(VEP.out.vcf)

    emit:
    vcf            = TABIX_BGZIPTABIX.out.tbi // channel: [ val(meta), vcf, tbi ]
    vep_reports = VEP.out.reports       //    path: *.html
    vep_version = VEP.out.version       //    path: *.version.txt
}
