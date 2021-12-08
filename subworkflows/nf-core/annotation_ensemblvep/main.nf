//
// Run VEP to annotate VCF files
//

include { ENSEMBLVEP                                } from '../../../modules/ensemblvep/main'
include { TABIX_BGZIPTABIX as ANNOTATION_BGZIPTABIX } from '../../../modules/tabix/bgziptabix/main'

workflow ANNOTATION_ENSEMBLVEP {
    take:
    vcf               // channel: [ val(meta), vcf ]
    vep_genome        //   value: which genome
    vep_species       //   value: which species
    vep_cache_version //   value: which cache version
    vep_cache         //    path: path_to_vep_cache (optionnal)

    main:
    ENSEMBLVEP(vcf, vep_genome, vep_species, vep_cache_version, vep_cache)
    ANNOTATION_BGZIPTABIX(ENSEMBLVEP.out.vcf)

    ch_versions = ENSEMBLVEP.out.versions.first().mix(ANNOTATION_BGZIPTABIX.out.versions.first())

    emit:
    vcf_tbi  = ANNOTATION_BGZIPTABIX.out.gz_tbi // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    reports  = ENSEMBLVEP.out.report            //    path: *.html
    versions = ch_versions                      //    path: versions.yml
}
