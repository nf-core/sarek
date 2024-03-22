//
// Validation against truth.
//

include { RTGTOOLS_VCFEVAL                 } from '../../modules/nf-core/rtgtools/vcfeval/main'

workflow VCF_VALIDATE_SMALL_VARIANTS {

    take:

    main:
        // input: of rtgtools/vcfeval
        // tuple val(meta), path(query_vcf), path(query_vcf_tbi), path(truth_vcf), path(truth_vcf_tbi), path(truth_bed), path(evaluation_bed)
        // tuple val(meta2), path(sdf)

        RTGTOOLS_VCFEVAL ( ch_vcfeval_in, ch_sdf )

    emit:
        versions        = ch_versions // channel: [ path(versions.yml) ]
}
