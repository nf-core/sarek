//
//
// PARABRICKS MUTECTCALLER: GPU-accelerated tumor-normal somatic variant calling
//

include { PARABRICKS_MUTECTCALLER } from '../../../modules/nf-core/parabricks/mutectcaller/main'
include { TABIX_TABIX              } from '../../../modules/nf-core/tabix/tabix/main'

workflow BAM_VARIANT_CALLING_SOMATIC_PARABRICKS_MUTECTCALLER {
    take:
    reads                  // channel: [ meta, normal_reads, normal_index, tumor_reads, tumor_index ] (BAM or CRAM)
    fasta                  // channel: [ meta, fasta ]
    fasta_fai              // channel: [ meta, fasta_fai ] - required for CRAM input
    panel_of_normals       // path: panel_of_normals or []
    panel_of_normals_tbi   // path: panel_of_normals_tbi or []
    intervals_bed_combined // channel: intervals or []

    main:
    versions = Channel.empty()

    // Rearrange to [ meta, tumor_reads, tumor_index, normal_reads, normal_index, intervals ]
    // Use single-param closure: when no_intervals, intervals_bed_combined is Channel.value([])
    // and combine passes the 5-element reads tuple as a single LinkedList item.
    ch_input = reads
        .combine(intervals_bed_combined)
        .map { row ->
            def meta         = row[0]
            def normal_reads = row[1]
            def normal_index = row[2]
            def tumor_reads  = row[3]
            def tumor_index  = row[4]
            def intervals    = row.size() > 5 ? row[5] : []
            [ meta, tumor_reads, tumor_index, normal_reads, normal_index, intervals ]
        }

    PARABRICKS_MUTECTCALLER(
        ch_input,
        fasta,
        fasta_fai,
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
