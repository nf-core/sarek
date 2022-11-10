//
// Run VEP to annotate VCF files
//

include { ENSEMBLVEP       } from '../../../modules/nf-core/ensemblvep/main'
include { TABIX_BGZIPTABIX } from '../../../modules/nf-core/tabix/bgziptabix/main'

workflow VCF_ANNOTATE_ENSEMBLVEP {
    take:
    vcf               // channel: [ val(meta), vcf ]
    fasta             //   value: fasta to use
    vep_genome        //   value: genome to use
    vep_species       //   value: species to use
    vep_cache_version //   value: cache version to use
    vep_cache         //    path: /path/to/vep/cache (optionnal)
    vep_extra_files   // channel: [ file1, file2...] (optionnal)

    main:
    ch_versions = Channel.empty()

    ENSEMBLVEP(vcf, vep_genome, vep_species, vep_cache_version, vep_cache, fasta, vep_extra_files)
    TABIX_BGZIPTABIX(ENSEMBLVEP.out.vcf)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(ENSEMBLVEP.out.versions.first())
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions.first())

    emit:
    vcf_tbi  = TABIX_BGZIPTABIX.out.gz_tbi // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    json     = ENSEMBLVEP.out.json
    tab      = ENSEMBLVEP.out.tab
    reports  = ENSEMBLVEP.out.report       //    path: *.html
    versions = ch_versions                 //    path: versions.yml
}
