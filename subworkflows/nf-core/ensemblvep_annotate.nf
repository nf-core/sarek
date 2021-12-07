//
// Run VEP to annotate VCF files
//

include { ENSEMBLVEP                              } from '../../modules/nf-core/modules/ensemblvep/main'
include { TABIX_BGZIPTABIX as BGZIPTABIX_ANNOTATE } from '../../modules/nf-core/modules/tabix/bgziptabix/main'

workflow ENSEMBLVEP_ANNOTATE {
    take:
    vcf               // channel: [ val(meta), vcf ]
    vep_genome        //   value: which genome
    vep_species       //   value: which species
    vep_cache_version //   value: which cache version
    vep_cache         //    path: path_to_vep_cache (optionnal)

    main:
    ENSEMBLVEP(vcf, vep_genome, vep_species, vep_cache_version, vep_cache)
    BGZIPTABIX_ANNOTATE(ENSEMBLVEP.out.vcf)

    emit:
    vcf_tbi  = BGZIPTABIX_ANNOTATE.out.gz_tbi // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    reports  = ENSEMBLVEP.out.report          //    path: *.html
    versions = ENSEMBLVEP.out.versions        //    path: versions.yml
}
