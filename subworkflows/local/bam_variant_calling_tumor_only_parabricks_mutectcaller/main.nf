//
//
// PARABRICKS MUTECTCALLER: GPU-accelerated tumor-only somatic variant calling
//

include { PARABRICKS_MUTECTCALLER               } from '../../../modules/nf-core/parabricks/mutectcaller/main'
include { HTSLIB_BGZIPTABIX as TABIX_BGZIPTABIX } from '../../../modules/nf-core/htslib/bgziptabix/main'

workflow BAM_VARIANT_CALLING_TUMOR_ONLY_PARABRICKS_MUTECTCALLER {
    take:
    reads                  // channel: [ meta, tumor_reads, tumor_index ] (BAM or CRAM)
    fasta                  // channel: [ meta, fasta ]
    fasta_fai              // channel: [ meta, fasta_fai ] - required for CRAM input
    panel_of_normals       // path: panel_of_normals or []
    panel_of_normals_tbi   // path: panel_of_normals_tbi or []
    intervals_bed_combined // channel: intervals or []

    main:
    // Rearrange to [ meta, tumor_reads, tumor_index, normal_reads, normal_index, intervals ]
    // Use single-param closure: when no_intervals, intervals_bed_combined is Channel.value([])
    // and combine passes the 3-element reads tuple as a single LinkedList item.
    ch_input = reads
        .combine(intervals_bed_combined)
        .map { row ->
            def meta        = row[0]
            def tumor_reads = row[1]
            def tumor_index = row[2]
            def intervals   = row.size() > 3 ? row[3] : []
            [ meta, tumor_reads, tumor_index, [], [], intervals ]
        }

    fasta_with_fai = fasta.combine(fasta_fai.map { _meta, fai -> fai }).collect()

    PARABRICKS_MUTECTCALLER(
        ch_input,
        fasta_with_fai,
        panel_of_normals,
        panel_of_normals_tbi,
    )


    TABIX_BGZIPTABIX(PARABRICKS_MUTECTCALLER.out.vcf.map { meta, vcf -> [ meta, vcf, [], [] ] }, 'compress', true, 'vcf')

    emit:
    vcf   = TABIX_BGZIPTABIX.out.output.map { meta, vcf -> [ meta + [ variantcaller: 'parabricks_mutectcaller' ], vcf ] }
    tbi   = TABIX_BGZIPTABIX.out.index.map  { meta, tbi -> [ meta + [ variantcaller: 'parabricks_mutectcaller' ], tbi ] }
    stats = PARABRICKS_MUTECTCALLER.out.stats
}
