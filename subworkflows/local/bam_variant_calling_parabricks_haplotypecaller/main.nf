//
// PARABRICKS HAPLOTYPECALLER germline variant calling (GPU-accelerated)
//

include { PARABRICKS_HAPLOTYPECALLER } from '../../../modules/nf-core/parabricks/haplotypecaller/main'
include { TABIX_BGZIPTABIX            } from '../../../modules/nf-core/tabix/bgziptabix/main'

workflow BAM_VARIANT_CALLING_PARABRICKS_HAPLOTYPECALLER {
    take:
    cram                    // channel: [mandatory] [ meta, cram, crai ]
    fasta                   // channel: [mandatory] [ meta, fasta ]
    intervals_bed_combined  // channel: [optional]  [] or [ intervals.bed ]

    main:
    // Combine each sample with the (optional) intervals list
    // intervals_bed_combined emits [] (no intervals) or [file] (one combined BED)
    // When no_intervals, the empty list contributes 0 elements to the combined tuple,
    // so use it.size() to safely extract the optional 4th element.
    cram_intervals = cram
        .combine(intervals_bed_combined)
        .map { it ->
            def (meta, cram_, crai) = it
            def intervals_ = it.size() > 3 ? it[3] : []
            [ meta + [ variantcaller:'parabricks_haplotypecaller' ], cram_, crai, intervals_ ]
        }

    PARABRICKS_HAPLOTYPECALLER(cram_intervals, fasta)

    // Compress and index the uncompressed VCF output
    TABIX_BGZIPTABIX(PARABRICKS_HAPLOTYPECALLER.out.vcf)

    vcf = TABIX_BGZIPTABIX.out.gz_index.map { meta, vcf_, _tbi -> [ meta, vcf_ ] }
    tbi = TABIX_BGZIPTABIX.out.gz_index.map { meta, _vcf, tbi  -> [ meta, tbi  ] }

    emit:
    vcf      // channel: [ val(meta), vcf.gz ]
    tbi      // channel: [ val(meta), vcf.gz.tbi ]
}
