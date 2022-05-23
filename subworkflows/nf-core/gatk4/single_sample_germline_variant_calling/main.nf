include { GATK4_CNNSCOREVARIANTS      as CNNSCOREVARIANTS      } from '../../../../modules/nf-core/modules/gatk4/cnnscorevariants/main'
include { GATK4_FILTERVARIANTTRANCHES as FILTERVARIANTTRANCHES } from  '../../../../modules/nf-core/modules/gatk4/filtervarianttranches/main'

workflow GATK_SINGLE_SAMPLE_GERMLINE_VARIANT_CALLING{

    take:
    vcf // meta, vcf, tbi
    fasta
    fasta_fai
    dict
    germline_resource //TODO more resources here, dbsnp, omni, 1000g, hapmap
    germline_resource_tbi

    main:

    ch_versions = Channel.empty()

    CNNSCOREVARIANTS(vcf.map{meta, vcf,tbi -> [meta,vcf,tbi,[],[]]}, fasta, fasta_fai, dict, [], [])
    FILTERVARIANTTRANCHES(CNNSCOREVARIANTS.out.vcf.join(CNNSCOREVARIANTS.out.tbi).map{meta, vcf,tbi -> [meta,vcf,tbi,[]]}, germline_resource, germline_resource_tbi, fasta, fasta_fai, dict)

    ch_versions = ch_versions.mix(CNNSCOREVARIANTS.out.versions)
    ch_versions = ch_versions.mix(FILTERVARIANTTRANCHES.out.versions)

    emit:
    versions = ch_versions
}
