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
    // Reshape intervals_bed_combined ([] or [file]) into a fixed-size 1-tuple ([] or file),
    // mirroring the deepvariant intervals channel, so combine() always yields a 4-tuple
    // and no destructuring/size-check is needed (there is no interval number here,
    // as parabricks processes everything in a single fat process).
    intervals_bed = intervals_bed_combined.map { bed -> [ bed ? bed[0] : [] ] }

    cram_intervals = cram.combine(intervals_bed)

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
