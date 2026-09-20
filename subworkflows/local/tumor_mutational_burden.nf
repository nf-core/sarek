include { BCFTOOLS_NORM } from '../../modules/nf-core/bcftools/norm/main'
include { TMB           } from '../../modules/local/tmb/main'

workflow TUMOR_MUTATIONAL_BURDEN {

    take:
    vcf
    fasta
    target_bed

    main:

    ch_versions = Channel.empty()

    BCFTOOLS_NORM(vcf, fasta)

     tmb_in = BCFTOOLS_NORM.out.vcf.map{ meta,vcf ->
         dbConfig = file("${projectDir}/${WorkflowSarek.getTMBdatabase(meta.annotation)}", checkIfExists: true)
         varConfig = file("${projectDir}/${WorkflowSarek.getTMBvariantcaller(meta.variantcaller)}", checkIfExists: true)

        [meta,vcf,dbConfig,varConfig]
    }

    TMB(tmb_in, target_bed)

    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)
    ch_versions = ch_versions.mix(TMB.out.versions)

    emit:
    versions                = ch_versions
}
