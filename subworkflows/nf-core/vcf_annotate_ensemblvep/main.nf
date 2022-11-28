//
// Run VEP to annotate VCF files
//

include { ENSEMBLVEP  } from '../../../modules/nf-core/ensemblvep/main'
include { TABIX_TABIX } from '../../../modules/nf-core/tabix/tabix/main'

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

    ENSEMBLVEP(vcf, vep_genome, vep_species, vep_cache_version, vep_cache, fasta, vep_extra_files)
    TABIX_TABIX(ENSEMBLVEP.out.vcf)

    ch_vcf_tbi = ENSEMBLVEP.out.vcf.join(TABIX_TABIX.out.tbi)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(ENSEMBLVEP.out.versions)
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    emit:
    vcf_tbi  = ch_vcf_tbi                  // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    json     = ENSEMBLVEP.out.json         // channel: [ val(meta), json ]
    tab      = ENSEMBLVEP.out.tab          // channel: [ val(meta), tab ]
    reports  = ENSEMBLVEP.out.report       //    path: *.html
    versions = ch_versions                 //    path: versions.yml
}
