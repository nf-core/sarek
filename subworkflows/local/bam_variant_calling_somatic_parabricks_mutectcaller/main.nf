//
// PARABRICKS MUTECTCALLER: GPU-accelerated tumor-normal somatic variant calling
//

include { HTSLIB_BGZIPTABIX as TABIX_PARABRICKS_MUTECTCALLER } from '../../../modules/nf-core/htslib/bgziptabix'
include { PARABRICKS_MUTECTCALLER                            } from '../../../modules/nf-core/parabricks/mutectcaller'

workflow BAM_VARIANT_CALLING_SOMATIC_PARABRICKS_MUTECTCALLER {
    take:
    reads                  // channel: [ meta, normal_reads, normal_index, tumor_reads, tumor_index ] (BAM or CRAM)
    fasta                  // channel: [ meta, fasta ]
    fasta_fai              // channel: [ meta, fasta_fai ] - required for CRAM input
    panel_of_normals       // channel: [ panel_of_normals ] or []
    panel_of_normals_tbi   // channel: [ panel_of_normals_tbi ] or []
    intervals_bed_combined // channel: [ intervals.bed ] or []

    main:
    // Rearrange to [ meta, tumor_reads, tumor_index, normal_reads, normal_index, intervals ]
    // Use single-param closure: when no_intervals, intervals_bed_combined is channel.value([])
    // and combine passes the 5-element reads tuple as a single LinkedList item.
    // Parabricks runs a single job per pair, so the combined intervals bed (or none) is used, no scatter/gather.
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

    fasta_with_fai = fasta.combine(fasta_fai.map { _meta, fai -> fai }).collect()

    PARABRICKS_MUTECTCALLER(
        ch_input,
        fasta_with_fai,
        panel_of_normals,
        panel_of_normals_tbi,
    )

    // The Parabricks VCF is already bgzip-compressed: only link it under its final name and index it (no re-compression)
    TABIX_PARABRICKS_MUTECTCALLER(PARABRICKS_MUTECTCALLER.out.vcf.map { meta, vcf_ -> [ meta, vcf_, [], [] ] }, 'compress', true, 'vcf')

    emit:
    vcf   = TABIX_PARABRICKS_MUTECTCALLER.out.output.map { meta, vcf_ -> [ meta + [ variantcaller: 'parabricks_mutectcaller' ], vcf_ ] } // channel: [ meta, vcf ]
    tbi   = TABIX_PARABRICKS_MUTECTCALLER.out.index.map  { meta, tbi_ -> [ meta + [ variantcaller: 'parabricks_mutectcaller' ], tbi_ ] } // channel: [ meta, tbi ]
    stats = PARABRICKS_MUTECTCALLER.out.stats                                                                                          // channel: [ meta, vcf.gz.stats ]
}
