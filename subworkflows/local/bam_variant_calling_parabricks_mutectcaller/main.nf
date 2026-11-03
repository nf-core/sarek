//
// PARABRICKS MUTECTCALLER: GPU-accelerated tumour-normal or tumour-only somatic variant calling, filtered with GATK on CPU
//

include { VCF_VARIANT_FILTERING_MUTECT2 } from '../vcf_variant_filtering_mutect2'
include { PARABRICKS_MUTECTCALLER       } from '../../../modules/nf-core/parabricks/mutectcaller'

workflow BAM_VARIANT_CALLING_PARABRICKS_MUTECTCALLER {
    take:
    cram                            // channel: [mandatory] [ meta, normal_cram, normal_crai, tumour_cram, tumour_crai ], normal_cram and normal_crai are [] for tumour-only
    fasta                           // channel: [mandatory] [ meta, fasta ]
    fasta_fai                       // channel: [mandatory] [ meta, fasta_fai ]
    dict                            // channel: [mandatory] [ meta, dict ]
    germline_resource               // channel: [optional] [ germline_resource ] or []
    germline_resource_tbi           // channel: [optional] [ germline_resource_tbi ] or []
    contamination_resource          // channel: [optional] [ contamination_resource ] or [], common biallelic SNPs for getpileupsummaries
    contamination_resource_tbi      // channel: [optional] [ contamination_resource_tbi ] or []
    panel_of_normals                // channel: [optional] [ panel_of_normals ] or []
    panel_of_normals_tbi            // channel: [optional] [ panel_of_normals_tbi ] or []
    intervals_bed_combined          // channel: [optional] [] or [ intervals.bed ]
    intervals                       // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals, scatter for getpileupsummaries
    contamination_resource_provided // boolean: [mandatory] whether a contamination resource was provided, so contamination tables are expected

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

    VCF_VARIANT_FILTERING_MUTECT2(
        calls,
        cram_filter,
        intervals,
        fasta,
        fasta_fai,
        dict,
        contamination_resource,
        contamination_resource_tbi,
        contamination_resource_provided,
        false,
    )

    vcf                 = VCF_VARIANT_FILTERING_MUTECT2.out.vcf
    tbi                 = VCF_VARIANT_FILTERING_MUTECT2.out.tbi
    stats               = PARABRICKS_MUTECTCALLER.out.stats
    stats_filtered      = VCF_VARIANT_FILTERING_MUTECT2.out.stats_filtered
    artifact_priors     = VCF_VARIANT_FILTERING_MUTECT2.out.artifact_priors
    contamination_table = VCF_VARIANT_FILTERING_MUTECT2.out.contamination_table
    segmentation_table  = VCF_VARIANT_FILTERING_MUTECT2.out.segmentation_table
    pileup_table_normal = VCF_VARIANT_FILTERING_MUTECT2.out.pileup_table_normal
    pileup_table_tumor  = VCF_VARIANT_FILTERING_MUTECT2.out.pileup_table_tumor

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
