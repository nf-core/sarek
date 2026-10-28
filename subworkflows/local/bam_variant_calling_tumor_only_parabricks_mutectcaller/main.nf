//
// PARABRICKS MUTECTCALLER: GPU-accelerated tumour-only somatic variant calling, filtered with GATK on CPU
//

include { GATK4_CALCULATECONTAMINATION as MUTECTCALLER_CALCULATECONTAMINATION } from '../../../modules/nf-core/gatk4/calculatecontamination'
include { GATK4_FILTERMUTECTCALLS as MUTECTCALLER_FILTERMUTECTCALLS           } from '../../../modules/nf-core/gatk4/filtermutectcalls'
include { GATK4_GETPILEUPSUMMARIES as MUTECTCALLER_GETPILEUPSUMMARIES         } from '../../../modules/nf-core/gatk4/getpileupsummaries'
include { PARABRICKS_MUTECTCALLER                                             } from '../../../modules/nf-core/parabricks/mutectcaller'

workflow BAM_VARIANT_CALLING_TUMOR_ONLY_PARABRICKS_MUTECTCALLER {
    take:
    cram                   // channel: [mandatory] [ meta, tumour_cram, tumour_crai ]
    fasta                  // channel: [mandatory] [ meta, fasta ]
    fasta_fai              // channel: [mandatory] [ meta, fasta_fai ]
    dict                   // channel: [mandatory] [ meta, dict ]
    germline_resource      // channel: [optional] [ germline_resource ] or []
    germline_resource_tbi  // channel: [optional] [ germline_resource_tbi ] or []
    panel_of_normals       // channel: [optional] [ panel_of_normals ] or []
    panel_of_normals_tbi   // channel: [optional] [ panel_of_normals_tbi ] or []
    intervals_bed_combined // channel: [optional] [] or [ intervals.bed ]

    main:
    intervals_bed = intervals_bed_combined.map { bed -> [bed ? bed[0] : []] }

    cram_intervals = cram
        .combine(intervals_bed)
        .map { meta, tumor_cram, tumor_crai, interval ->
            [meta, tumor_cram, tumor_crai, [], [], interval]
        }

    fasta_with_fai = fasta.combine(fasta_fai.map { _meta, fai -> fai }).collect()

    PARABRICKS_MUTECTCALLER(cram_intervals, fasta_with_fai, panel_of_normals, panel_of_normals_tbi)

    // Without a germline resource the pileup and contamination steps are skipped
    germline_resource_pileup     = germline_resource.filter { resource -> resource != [] }
    germline_resource_pileup_tbi = germline_resource_tbi.filter { resource_tbi -> resource_tbi != [] }

    MUTECTCALLER_GETPILEUPSUMMARIES(
        cram.combine(intervals_bed),
        fasta,
        fasta_fai,
        dict,
        germline_resource_pileup,
        germline_resource_pileup_tbi,
    )

    MUTECTCALLER_CALCULATECONTAMINATION(MUTECTCALLER_GETPILEUPSUMMARIES.out.table.map { meta, table -> [meta, table, []] })

    // Contamination and segmentation tables only exist with a germline resource; no orientation priors are available
    vcf_to_filter = PARABRICKS_MUTECTCALLER.out.vcf
        .join(PARABRICKS_MUTECTCALLER.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .join(PARABRICKS_MUTECTCALLER.out.stats, failOnDuplicate: true, failOnMismatch: true)
        .join(MUTECTCALLER_CALCULATECONTAMINATION.out.segmentation, failOnDuplicate: true, remainder: true)
        .join(MUTECTCALLER_CALCULATECONTAMINATION.out.contamination, failOnDuplicate: true, remainder: true)
        .filter { row -> row[1] != null }
        .map { meta, vcf, tbi, stats, segmentation, contamination -> [meta, vcf, tbi, stats, [], segmentation ?: [], contamination ?: [], []] }

    MUTECTCALLER_FILTERMUTECTCALLS(vcf_to_filter, fasta, fasta_fai, dict)

    vcf = MUTECTCALLER_FILTERMUTECTCALLS.out.vcf.map { meta, vcf -> [meta + [variantcaller: 'parabricks_mutectcaller'], vcf] }

    tbi = MUTECTCALLER_FILTERMUTECTCALLS.out.tbi.map { meta, tbi -> [meta + [variantcaller: 'parabricks_mutectcaller'], tbi] }

    stats               = PARABRICKS_MUTECTCALLER.out.stats
    stats_filtered      = MUTECTCALLER_FILTERMUTECTCALLS.out.stats
    contamination_table = MUTECTCALLER_CALCULATECONTAMINATION.out.contamination
    segmentation_table  = MUTECTCALLER_CALCULATECONTAMINATION.out.segmentation
    pileup_table        = MUTECTCALLER_GETPILEUPSUMMARIES.out.table

    emit:
    vcf                 // channel: [ meta, vcf.gz ]
    tbi                 // channel: [ meta, vcf.gz.tbi ]
    stats               // channel: [ meta, vcf.gz.stats ]
    stats_filtered      // channel: [ meta, filteringStats.tsv ]
    contamination_table // channel: [ meta, contamination.table ]
    segmentation_table  // channel: [ meta, segmentation.table ]
    pileup_table        // channel: [ meta, pileups.table ]
}
