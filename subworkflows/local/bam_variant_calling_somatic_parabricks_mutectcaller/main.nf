//
// PARABRICKS MUTECTCALLER: GPU-accelerated tumour-normal somatic variant calling
//

include { PARABRICKS_MUTECTCALLER } from '../../../modules/nf-core/parabricks/mutectcaller'

workflow BAM_VARIANT_CALLING_SOMATIC_PARABRICKS_MUTECTCALLER {
    take:
    cram                   // channel: [mandatory] [ meta, normal_cram, normal_crai, tumour_cram, tumour_crai ]
    fasta                  // channel: [mandatory] [ meta, fasta ]
    fasta_fai              // channel: [mandatory] [ meta, fasta_fai ]
    panel_of_normals       // channel: [optional] [ panel_of_normals ] or []
    panel_of_normals_tbi   // channel: [optional] [ panel_of_normals_tbi ] or []
    intervals_bed_combined // channel: [optional] [] or [ intervals.bed ]

    main:
    intervals_bed = intervals_bed_combined.map { bed -> [bed ? bed[0] : []] }

    cram_intervals = cram
        .combine(intervals_bed)
        .map { meta, normal_cram, normal_crai, tumor_cram, tumor_crai, interval ->
            [meta, tumor_cram, tumor_crai, normal_cram, normal_crai, interval]
        }

    fasta_with_fai = fasta.combine(fasta_fai.map { _meta, fai -> fai }).collect()

    PARABRICKS_MUTECTCALLER(cram_intervals, fasta_with_fai, panel_of_normals, panel_of_normals_tbi)

    vcf = PARABRICKS_MUTECTCALLER.out.vcf.map { meta, vcf -> [meta + [variantcaller: 'parabricks_mutectcaller'], vcf] }

    tbi = PARABRICKS_MUTECTCALLER.out.tbi.map { meta, tbi -> [meta + [variantcaller: 'parabricks_mutectcaller'], tbi] }

    stats = PARABRICKS_MUTECTCALLER.out.stats

    emit:
    vcf   // channel: [ meta, vcf.gz ]
    tbi   // channel: [ meta, vcf.gz.tbi ]
    stats // channel: [ meta, vcf.gz.stats ]
}
