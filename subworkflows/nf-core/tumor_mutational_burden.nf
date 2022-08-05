include { BCFTOOLS_NORM } from '../../modules/nf-core/modules/bcftools/norm/main'
include { TMB           } from '../../modules/nf-core/modules/tmb/main'

workflow TUMOR_MUTATIONAL_BURDEN {

    take:
    vcf
    fasta
    target_bed

    main:

    ch_versions = Channel.empty()

    vcf.view()
    //bcftools norm -f FASTA -m- -o file_norm.vcf file
    BCFTOOLS_NORM(vcf, fasta)

    BCFTOOLS_NORM.out.vcf.view()

    // tmb_in = BCFTOOLS_NORM.out.vcf.map{ meta,vcf ->
        // dbConfig = file(WorkflowSarek.getTMBdatabase(meta.annotation), checkIfExists: true)
        // varConfig = file(WorkflowSarek.getTMBvariantcaller(meta.variantcaller), checkIfExists: true)
//
        // [meta,vcf,dbConfig,varConfig]
    // }

    //TMB(tmb_in, target_bed)

    emit:
    versions                = ch_versions
}
