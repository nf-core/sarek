//
// PARABRICKS MUTECTCALLER: GPU-accelerated tumour-normal or tumour-only somatic variant calling, filtered with GATK on CPU
//

include { BAM_VCF_FILTER_MUTECT_CALLS_GATK } from '../bam_vcf_filter_mutect_calls_gatk/main'
include { PARABRICKS_MUTECTCALLER          } from '../../../modules/nf-core/parabricks/mutectcaller'

workflow BAM_VARIANT_CALLING_PARABRICKS_MUTECTCALLER {
    take:
    cram                       // channel: [mandatory] [ meta, normal_cram, normal_crai, tumour_cram, tumour_crai ], normal_cram and normal_crai are [] for tumour-only
    fasta                      // channel: [mandatory] [ meta, fasta ]
    fasta_fai                  // channel: [mandatory] [ meta, fasta_fai ]
    dict                       // channel: [mandatory] [ meta, dict ]
    germline_resource          // channel: [optional] [ germline_resource ] or []
    germline_resource_tbi      // channel: [optional] [ germline_resource_tbi ] or []
    panel_of_normals           // channel: [optional] [ panel_of_normals ] or []
    panel_of_normals_tbi       // channel: [optional] [ panel_of_normals_tbi ] or []
    intervals_bed_combined     // channel: [optional] [] or [ intervals.bed ]
    intervals                  // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals, scatter for getpileupsummaries
    germline_resource_provided // boolean: [mandatory] whether a germline resource was provided, so contamination tables are expected

    main:
    intervals_bed = intervals_bed_combined.map { bed -> [bed ? bed[0] : []] }

    cram_intervals = cram
        .combine(intervals_bed)
        .map { meta, normal_cram, normal_crai, tumor_cram, tumor_crai, interval ->
            [meta, tumor_cram, tumor_crai, normal_cram, normal_crai, interval]
        }

    fasta_with_fai = fasta.combine(fasta_fai.map { _meta, fai -> fai }).collect()

    PARABRICKS_MUTECTCALLER(cram_intervals, fasta_with_fai, panel_of_normals, panel_of_normals_tbi, germline_resource, germline_resource_tbi)

    calls = PARABRICKS_MUTECTCALLER.out.vcf
        .join(PARABRICKS_MUTECTCALLER.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .join(PARABRICKS_MUTECTCALLER.out.stats, failOnDuplicate: true, failOnMismatch: true)
        .join(PARABRICKS_MUTECTCALLER.out.f1r2, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, vcf_, tbi_, stats_, f1r2 -> [meta + [variantcaller: 'parabricks_mutectcaller'], vcf_, tbi_, stats_, f1r2] }

    cram_filter = cram.map { meta, normal_cram, normal_crai, tumor_cram, tumor_crai ->
        [meta + [variantcaller: 'parabricks_mutectcaller'], normal_cram, normal_crai, tumor_cram, tumor_crai]
    }

    BAM_VCF_FILTER_MUTECT_CALLS_GATK(
        calls,
        cram_filter,
        intervals,
        fasta,
        fasta_fai,
        dict,
        germline_resource,
        germline_resource_tbi,
        germline_resource_provided,
        false,
    )

    vcf                 = BAM_VCF_FILTER_MUTECT_CALLS_GATK.out.vcf
    tbi                 = BAM_VCF_FILTER_MUTECT_CALLS_GATK.out.tbi
    stats               = PARABRICKS_MUTECTCALLER.out.stats
    stats_filtered      = BAM_VCF_FILTER_MUTECT_CALLS_GATK.out.stats_filtered
    artifact_priors     = BAM_VCF_FILTER_MUTECT_CALLS_GATK.out.artifact_priors
    contamination_table = BAM_VCF_FILTER_MUTECT_CALLS_GATK.out.contamination_table
    segmentation_table  = BAM_VCF_FILTER_MUTECT_CALLS_GATK.out.segmentation_table
    pileup_table_normal = BAM_VCF_FILTER_MUTECT_CALLS_GATK.out.pileup_table_normal
    pileup_table_tumor  = BAM_VCF_FILTER_MUTECT_CALLS_GATK.out.pileup_table_tumor

    emit:
    vcf                 // channel: [ meta, vcf.gz ]
    tbi                 // channel: [ meta, vcf.gz.tbi ]
    stats               // channel: [ meta, vcf.gz.stats ]
    stats_filtered      // channel: [ meta, filteringStats.tsv ]
    artifact_priors     // channel: [ meta, artifactprior.tar.gz ]
    contamination_table // channel: [ meta, contamination.table ]
    segmentation_table  // channel: [ meta, segmentation.table ]
    pileup_table_normal // channel: [ meta, pileups.table ]
    pileup_table_tumor  // channel: [ meta, pileups.table ]
}
