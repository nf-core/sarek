//
// PARABRICKS DEEPVARIANT germline variant calling (GPU-accelerated)
//

include { PARABRICKS_DEEPVARIANT                              } from '../../../modules/nf-core/parabricks/deepvariant/main'
include { HTSLIB_BGZIPTABIX as TABIX_VC_PARABRICKS_DEEPVARIANT } from '../../../modules/nf-core/htslib/bgziptabix/main'

workflow BAM_VARIANT_CALLING_PARABRICKS_DEEPVARIANT {
    take:
    cram                    // channel: [mandatory] [ meta, cram, crai ]
    fasta                   // channel: [mandatory] [ meta, fasta ]
    intervals_bed_combined  // channel: [optional]  [] or [ intervals.bed ]

    main:
    // Combine each sample with the (optional) intervals list
    // intervals_bed_combined emits [] (no intervals) or [file] (one combined BED)
    // When no_intervals, the empty list contributes 0 elements to the combined tuple,
    // so check the tuple size to safely extract the optional 4th element.
    cram_intervals = cram
        .combine(intervals_bed_combined)
        .map { cram_combined ->
            def (meta, cram_, crai) = cram_combined
            def intervals_ = cram_combined.size() > 3 ? cram_combined[3] : []
            [ meta, cram_, crai, intervals_ ]
        }

    PARABRICKS_DEEPVARIANT(
        cram_intervals,
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
