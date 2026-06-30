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
    versions = Channel.empty()

    // Combine each sample with the (optional) intervals list
    // intervals_bed_combined emits [] (no intervals) or [file] (one combined BED)
    cram_intervals = cram
        .combine(intervals_bed_combined)
        .map { meta, cram_, crai, intervals_ ->
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
    versions // channel: versions.yml
}
