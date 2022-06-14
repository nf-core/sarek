include { GATK4_CNNSCOREVARIANTS      as CNNSCOREVARIANTS      } from '../../../../modules/nf-core/modules/gatk4/cnnscorevariants/main'
include { GATK4_FILTERVARIANTTRANCHES as FILTERVARIANTTRANCHES } from '../../../../modules/nf-core/modules/gatk4/filtervarianttranches/main'
include { GATK4_HAPLOTYPECALLER       as HAPLOTYPECALLER       } from '../../../../modules/nf-core/modules/gatk4/haplotypecaller/main'

workflow GATK_SINGLE_SAMPLE_GERMLINE_VARIANT_CALLING{

    take:
    cram // meta, cram, crai, intervals
    fasta
    fasta_fai
    dict
    known_sites
    known_sites_tbi
    dbsnp
    dbsnp_tbi

    main:

    ch_versions = Channel.empty()

    HAPLOTYPECALLER(
        cram,
        fasta,
        fasta_fai,
        dict,
        dbsnp,
        dbsnp_tbi)

    //Adding intervals
    CNNSCOREVARIANTS(HAPLOTYPECALLER.out.map{meta, vcf,tbi -> [meta,vcf,tbi,[],[]]}, fasta, fasta_fai, dict, [], [])

    //Adding intervals
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
    vcf = FILTERVARIANTTRANCHES.out.vcf
}
