include { GATK4_CNNSCOREVARIANTS      } from '../../../../modules/nf-core/modules/gatk4/cnnscorevariants/main'
include { GATK4_FILTERVARIANTTRANCHES } from  '../../../../modules/nf-core/modules/gatk4/filtervarianttranches/main'

workflow GATK_SINGLE_SAMPLE_GERMLINE_VARIANT_CALLING{

    take:
    vcf

    main:

    versions = Channel.empty()

    emit:
    versions
}
