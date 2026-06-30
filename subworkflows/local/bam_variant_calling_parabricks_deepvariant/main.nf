//
// PARABRICKS DEEPVARIANT germline variant calling (GPU-accelerated)
//

include { PARABRICKS_DEEPVARIANT                          } from '../../../modules/nf-core/parabricks/deepvariant/main'
include { TABIX_TABIX as TABIX_VC_PARABRICKS_DEEPVARIANT  } from '../../../modules/nf-core/tabix/tabix'

workflow BAM_VARIANT_CALLING_PARABRICKS_DEEPVARIANT {
    take:
    bam       // channel: [mandatory] [ meta, bam, bai ]
    fasta     // channel: [mandatory] [ meta, fasta ]

    main:
    versions = Channel.empty()

    PARABRICKS_DEEPVARIANT(
        bam.map { meta, bam_, bai -> [ meta, bam_, bai, [] ] },
        fasta
    )

    TABIX_VC_PARABRICKS_DEEPVARIANT(PARABRICKS_DEEPVARIANT.out.vcf)

    vcf = PARABRICKS_DEEPVARIANT.out.vcf
        .map { meta, vcf_ -> [ meta + [ variantcaller:'parabricks_deepvariant' ], vcf_ ] }

    tbi = TABIX_VC_PARABRICKS_DEEPVARIANT.out.tbi
        .map { meta, tbi_ -> [ meta + [ variantcaller:'parabricks_deepvariant' ], tbi_ ] }

    versions = versions.mix(TABIX_VC_PARABRICKS_DEEPVARIANT.out.versions)

    emit:
    vcf      // channel: [ meta, vcf.gz ]
    tbi      // channel: [ meta, vcf.gz.tbi ]
    versions // channel: [ versions.yml ]
}
