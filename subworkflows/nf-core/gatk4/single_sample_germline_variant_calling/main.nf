include { GATK4_CNNSCOREVARIANTS      as CNNSCOREVARIANTS      } from '../../../../modules/nf-core/modules/gatk4/cnnscorevariants/main'
include { GATK4_FILTERVARIANTTRANCHES as FILTERVARIANTTRANCHES } from  '../../../../modules/nf-core/modules/gatk4/filtervarianttranches/main'

workflow GATK_SINGLE_SAMPLE_GERMLINE_VARIANT_CALLING{

    take:
    vcf // meta, vcf, tbi
    fasta
    fasta_fai
    dict
       //TODO more resources here, dbsnp (done), omni, 1000g(done), hapmap
    known_sites
    known_sites_tbi

    main:

    ch_versions = Channel.empty()

    // intervals and merging?
    CNNSCOREVARIANTS(vcf.map{meta, vcf,tbi -> [meta,vcf,tbi,[],[]]}, fasta, fasta_fai, dict, [], [])

    FILTERVARIANTTRANCHES(CNNSCOREVARIANTS.out.vcf.join(CNNSCOREVARIANTS.out.tbi).map{meta, vcf,tbi -> [meta,vcf,tbi,[]]},
                            known_sites,
                            known_sites_tbi,
                            fasta,
                            fasta_fai,
                            dict)

    ch_versions = ch_versions.mix(CNNSCOREVARIANTS.out.versions)
    ch_versions = ch_versions.mix(FILTERVARIANTTRANCHES.out.versions)

    emit:
    versions = ch_versions
}
