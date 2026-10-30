//
// PARABRICKS MUTECTCALLER: GPU-accelerated tumour-normal or tumour-only somatic variant calling, filtered with GATK on CPU
//

include { GATK4_CALCULATECONTAMINATION as MUTECTCALLER_CALCULATECONTAMINATION       } from '../../../modules/nf-core/gatk4/calculatecontamination'
include { GATK4_FILTERMUTECTCALLS as MUTECTCALLER_FILTERMUTECTCALLS                 } from '../../../modules/nf-core/gatk4/filtermutectcalls'
include { GATK4_GETPILEUPSUMMARIES as MUTECTCALLER_GETPILEUPSUMMARIES_NORMAL        } from '../../../modules/nf-core/gatk4/getpileupsummaries'
include { GATK4_GETPILEUPSUMMARIES as MUTECTCALLER_GETPILEUPSUMMARIES_TUMOR         } from '../../../modules/nf-core/gatk4/getpileupsummaries'
include { GATK4_LEARNREADORIENTATIONMODEL as MUTECTCALLER_LEARNREADORIENTATIONMODEL } from '../../../modules/nf-core/gatk4/learnreadorientationmodel'
include { PARABRICKS_MUTECTCALLER                                                   } from '../../../modules/nf-core/parabricks/mutectcaller'

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

    MUTECTCALLER_LEARNREADORIENTATIONMODEL(PARABRICKS_MUTECTCALLER.out.f1r2)

    // Without a germline resource the pileup and contamination steps are skipped
    germline_resource_pileup     = germline_resource.filter { resource -> resource != [] }
    germline_resource_pileup_tbi = germline_resource_tbi.filter { resource_tbi -> resource_tbi != [] }

    // dropping the tumour-specific keys (tumor_id, contamination) lets a normal shared by several tumours be summarised once
    pileup = cram_intervals.multiMap { meta, tumor_cram, tumor_crai, normal_cram, normal_crai, interval ->
        tumor: [meta + [id: meta.tumor_id ?: meta.id], tumor_cram, tumor_crai, interval]
        normal: [meta - meta.subMap('tumor_id', 'contamination') + [id: meta.normal_id], normal_cram, normal_crai, interval]
    }

    pileup_normal = pileup.normal.filter { meta, _normal_cram, _normal_crai, _interval -> meta.normal_id }.unique()

    MUTECTCALLER_GETPILEUPSUMMARIES_NORMAL(pileup_normal, fasta, fasta_fai, dict, germline_resource_pileup, germline_resource_pileup_tbi)
    MUTECTCALLER_GETPILEUPSUMMARIES_TUMOR(pileup.tumor, fasta, fasta_fai, dict, germline_resource_pileup, germline_resource_pileup_tbi)

    pileup_table_tumor_branch = MUTECTCALLER_GETPILEUPSUMMARIES_TUMOR.out.table.branch { meta, _table ->
        paired: meta.normal_id
        tumor_only: true
    }

    pileup_normal_by_normal = MUTECTCALLER_GETPILEUPSUMMARIES_NORMAL.out.table.map { meta, table ->
        [[meta.patient, meta.normal_id], table]
    }

    calculatecontamination_in = pileup_table_tumor_branch.paired
        .map { meta, table ->
            [[meta.patient, meta.normal_id], meta + [id: "${meta.tumor_id}_vs_${meta.normal_id}".toString()], table]
        }
        .combine(pileup_normal_by_normal, by: 0)
        .map { _key, meta, tumor_table, normal_table -> [meta, tumor_table, normal_table] }
        .mix(pileup_table_tumor_branch.tumor_only.map { meta, table -> [meta, table, []] })

    MUTECTCALLER_CALCULATECONTAMINATION(calculatecontamination_in)

    vcf_with_priors = PARABRICKS_MUTECTCALLER.out.vcf
        .join(PARABRICKS_MUTECTCALLER.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .join(PARABRICKS_MUTECTCALLER.out.stats, failOnDuplicate: true, failOnMismatch: true)
        .join(MUTECTCALLER_LEARNREADORIENTATIONMODEL.out.artifactprior, failOnDuplicate: true, failOnMismatch: true)

    // Contamination and segmentation tables only exist with a germline resource; branching on the
    // boolean avoids a remainder join, which would hold every sample until the slowest caller finishes
    if (germline_resource_provided) {
        vcf_to_filter = vcf_with_priors
            .join(MUTECTCALLER_CALCULATECONTAMINATION.out.segmentation, failOnDuplicate: true, failOnMismatch: true)
            .join(MUTECTCALLER_CALCULATECONTAMINATION.out.contamination, failOnDuplicate: true, failOnMismatch: true)
            .map { meta, vcf, tbi, stats, artifactprior, segmentation, contamination ->
                [meta, vcf, tbi, stats, artifactprior, segmentation, contamination, []]
            }
    }
    else {
        vcf_to_filter = vcf_with_priors.map { meta, vcf, tbi, stats, artifactprior ->
            [meta, vcf, tbi, stats, artifactprior, [], [], []]
        }
    }

    MUTECTCALLER_FILTERMUTECTCALLS(vcf_to_filter, fasta, fasta_fai, dict)

    vcf = MUTECTCALLER_FILTERMUTECTCALLS.out.vcf.map { meta, vcf -> [meta + [variantcaller: 'parabricks_mutectcaller'], vcf] }

    tbi = MUTECTCALLER_FILTERMUTECTCALLS.out.tbi.map { meta, tbi -> [meta + [variantcaller: 'parabricks_mutectcaller'], tbi] }

    stats               = PARABRICKS_MUTECTCALLER.out.stats
    stats_filtered      = MUTECTCALLER_FILTERMUTECTCALLS.out.stats
    artifact_priors     = MUTECTCALLER_LEARNREADORIENTATIONMODEL.out.artifactprior
    contamination_table = MUTECTCALLER_CALCULATECONTAMINATION.out.contamination
    segmentation_table  = MUTECTCALLER_CALCULATECONTAMINATION.out.segmentation
    pileup_table_normal = MUTECTCALLER_GETPILEUPSUMMARIES_NORMAL.out.table
    pileup_table_tumor  = MUTECTCALLER_GETPILEUPSUMMARIES_TUMOR.out.table

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
