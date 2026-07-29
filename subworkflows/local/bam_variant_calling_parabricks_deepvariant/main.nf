//
// PARABRICKS DEEPVARIANT germline variant calling (GPU-accelerated)
//

include { PARABRICKS_DEEPVARIANT                              } from '../../../modules/nf-core/parabricks/deepvariant/main'
include { HTSLIB_BGZIPTABIX as TABIX_VC_PARABRICKS_DEEPVARIANT } from '../../../modules/nf-core/htslib/bgziptabix/main'

workflow BAM_VARIANT_CALLING_PARABRICKS_DEEPVARIANT {
    take:
    bam       // channel: [mandatory] [ meta, bam, bai ]
    fasta     // channel: [mandatory] [ meta, fasta ]

    main:
    PARABRICKS_DEEPVARIANT(
        bam.map { meta, bam_, bai -> [ meta, bam_, bai, [] ] },
        fasta
    )

    // Index the bgzip-compressed VCF output
    TABIX_VC_PARABRICKS_DEEPVARIANT(PARABRICKS_DEEPVARIANT.out.vcf.map { meta, vcf_ -> [ meta, vcf_, [], [] ] }, 'compress', true, 'vcf')

    vcf_tbi = TABIX_VC_PARABRICKS_DEEPVARIANT.out.output.join(TABIX_VC_PARABRICKS_DEEPVARIANT.out.index)

    vcf = vcf_tbi.map { meta, vcf_, _tbi -> [ meta + [ variantcaller:'parabricks_deepvariant' ], vcf_ ] }
    tbi = vcf_tbi.map { meta, _vcf, tbi_ -> [ meta + [ variantcaller:'parabricks_deepvariant' ], tbi_ ] }

    emit:
    vcf      // channel: [ meta, vcf.gz ]
    tbi      // channel: [ meta, vcf.gz.tbi ]
}
