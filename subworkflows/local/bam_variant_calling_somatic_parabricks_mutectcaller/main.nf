//
//
// PARABRICKS MUTECTCALLER: GPU-accelerated tumor-normal somatic variant calling
//

include { PARABRICKS_MUTECTCALLER } from '../../../modules/nf-core/parabricks/mutectcaller/main'
include { TABIX_TABIX              } from '../../../modules/nf-core/tabix/tabix/main'

workflow BAM_VARIANT_CALLING_SOMATIC_PARABRICKS_MUTECTCALLER {
    take:
    bam                    // channel: [ meta, normal_bam, normal_bai, tumor_bam, tumor_bai ]
    fasta                  // channel: [ meta, fasta ]
    panel_of_normals       // path: panel_of_normals or []
    panel_of_normals_tbi   // path: panel_of_normals_tbi or []
    intervals_bed_combined // channel: intervals or []

    main:
    versions = Channel.empty()

    // Rearrange to [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai, intervals ]
    ch_input = bam
        .combine(intervals_bed_combined)
        .map { meta, normal_bam, normal_bai, tumor_bam, tumor_bai, intervals ->
            [ meta, tumor_bam, tumor_bai, normal_bam, normal_bai, intervals ]
        }

    PARABRICKS_MUTECTCALLER(
        ch_input,
        fasta,
        panel_of_normals,
        panel_of_normals_tbi,
    )
    // PARABRICKS_MUTECTCALLER uses topic: versions — no out.versions to mix

    TABIX_TABIX(PARABRICKS_MUTECTCALLER.out.vcf)
    versions = versions.mix(TABIX_TABIX.out.versions)

    emit:
    vcf      = PARABRICKS_MUTECTCALLER.out.vcf.map { meta, vcf -> [ meta + [ variantcaller: 'parabricks_mutectcaller' ], vcf ] }
    tbi      = TABIX_TABIX.out.tbi.map           { meta, tbi -> [ meta + [ variantcaller: 'parabricks_mutectcaller' ], tbi ] }
    stats    = PARABRICKS_MUTECTCALLER.out.stats
    versions // channel: [ versions.yml ]
}
