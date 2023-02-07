//
// Run VEP to annotate VCF files
//

include { ENSEMBLVEP_VEP } from '../../../modules/nf-core/ensemblvep/vep/main'
include { TABIX_TABIX    } from '../../../modules/nf-core/tabix/tabix/main'

workflow VCF_ANNOTATE_ENSEMBLVEP {
    take:
    vcf               // channel: [ val(meta), vcf ]
    fasta             //   value: fasta to use (optionnal)
    vep_genome        //   value: genome to use
    vep_species       //   value: species to use
    vep_cache_version //   value: cache version to use
    vep_cache         //    path: /path/to/vep/cache (optionnal)
    vep_extra_files   // channel: [ file1, file2...] (optionnal)

    main:
    ch_versions = Channel.empty()

    ENSEMBLVEP_VEP(vcf, vep_genome, vep_species, vep_cache_version, vep_cache, fasta, vep_extra_files)
    TABIX_TABIX(ENSEMBLVEP_VEP.out.vcf)

    ch_vcf_tbi = ENSEMBLVEP_VEP.out.vcf.join(TABIX_TABIX.out.tbi)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions)
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    emit:
    vcf_tbi  = ch_vcf_tbi                  // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    json     = ENSEMBLVEP_VEP.out.json     // channel: [ val(meta), json ]
    tab      = ENSEMBLVEP_VEP.out.tab      // channel: [ val(meta), tab ]
    reports  = ENSEMBLVEP_VEP.out.report   // channel: [ *.html ]
    versions = ch_versions                 // channel: [ versions.yml ]
}
