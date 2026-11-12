//
// PARABRICKS HAPLOTYPECALLER germline variant calling (GPU-accelerated)
//

include { PARABRICKS_DBSNP           } from '../../../modules/nf-core/parabricks/dbsnp'
include { PARABRICKS_HAPLOTYPECALLER } from '../../../modules/nf-core/parabricks/haplotypecaller'

workflow BAM_VARIANT_CALLING_PARABRICKS_HAPLOTYPECALLER {
    take:
    cram                    // channel: [mandatory] [ meta, cram, crai ]
    fasta                   // channel: [mandatory] [ meta, fasta ]
    dbsnp                   // channel: [optional]  [] or [ dbsnp.vcf.gz ]
    dbsnp_tbi               // channel: [optional]  [] or [ dbsnp.vcf.gz.tbi ]
    intervals               // channel: [mandatory] [ interval.bed, num_intervals ], used to scatter joint genotyping
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
            [ meta + [ variantcaller:'parabricks_haplotypecaller' ], cram_, crai, intervals_ ]
        }

    PARABRICKS_HAPLOTYPECALLER(cram_intervals, fasta)

    // pbrun haplotypecaller has no --dbsnp: with --dbsnp it writes a plain VCF (ext.suffix), which pbrun dbsnp annotates, bgzips and indexes
    PARABRICKS_DBSNP(PARABRICKS_HAPLOTYPECALLER.out.vcf.combine(dbsnp).combine(dbsnp_tbi))

    vcf      = PARABRICKS_DBSNP.out.vcf.mix(PARABRICKS_HAPLOTYPECALLER.out.vcf_gz)
    tbi      = PARABRICKS_DBSNP.out.tbi.mix(PARABRICKS_HAPLOTYPECALLER.out.vcf_tbi)
    gvcf     = PARABRICKS_HAPLOTYPECALLER.out.gvcf
    gvcf_tbi = PARABRICKS_HAPLOTYPECALLER.out.gvcf_tbi

    // The gVCF covers all intervals in one file; GenomicsDBImport restricts each import to its own interval
    gvcf_tbi_intervals = gvcf.join(gvcf_tbi, failOnDuplicate: true, failOnMismatch: true)
        .combine(intervals)
        .map { meta, gvcf_, gvcf_tbi_, intervals_, num_intervals -> [ meta + [ num_intervals:num_intervals ], gvcf_, gvcf_tbi_, intervals_ ] }

    emit:
    vcf                 // channel: [ val(meta), vcf.gz ]
    tbi                 // channel: [ val(meta), vcf.gz.tbi ]
    gvcf                // channel: [ val(meta), g.vcf.gz ]
    gvcf_tbi            // channel: [ val(meta), g.vcf.gz.tbi ]
    gvcf_tbi_intervals  // channel: [ val(meta), g.vcf.gz, g.vcf.gz.tbi, interval.bed ], for joint genotyping
}
