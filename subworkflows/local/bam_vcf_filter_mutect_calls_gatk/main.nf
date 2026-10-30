//
// FILTER MUTECT CALLS: learnreadorientationmodel, getpileupsummaries, calculatecontamination and filtermutectcalls
// on the calls of GATK Mutect2 or Parabricks mutectcaller
//

include { GATK4_CALCULATECONTAMINATION as CALCULATECONTAMINATION       } from '../../../modules/nf-core/gatk4/calculatecontamination'
include { GATK4_FILTERMUTECTCALLS as FILTERMUTECTCALLS                 } from '../../../modules/nf-core/gatk4/filtermutectcalls'
include { GATK4_GATHERPILEUPSUMMARIES as GATHERPILEUPSUMMARIES_NORMAL  } from '../../../modules/nf-core/gatk4/gatherpileupsummaries'
include { GATK4_GATHERPILEUPSUMMARIES as GATHERPILEUPSUMMARIES_TUMOR   } from '../../../modules/nf-core/gatk4/gatherpileupsummaries'
include { GATK4_GETPILEUPSUMMARIES as GETPILEUPSUMMARIES_NORMAL        } from '../../../modules/nf-core/gatk4/getpileupsummaries'
include { GATK4_GETPILEUPSUMMARIES as GETPILEUPSUMMARIES_TUMOR         } from '../../../modules/nf-core/gatk4/getpileupsummaries'
include { GATK4_LEARNREADORIENTATIONMODEL as LEARNREADORIENTATIONMODEL } from '../../../modules/nf-core/gatk4/learnreadorientationmodel'

workflow BAM_VCF_FILTER_MUTECT_CALLS_GATK {
    take:
    calls                      // channel: [mandatory] [ meta, vcf, tbi, stats, f1r2 ] one row per pair, tumour-only sample or patient (joint), meta.variantcaller set
    cram                       // channel: [mandatory] [ meta, normal_cram, normal_crai, tumour_cram, tumour_crai ] one row per pair or tumour-only sample ([] for the normal)
    intervals                  // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals
    fasta                      // channel: [mandatory] [ meta, fasta ]
    fai                        // channel: [mandatory] [ meta, fai ]
    dict                       // channel: [mandatory] [ meta, dict ]
    germline_resource          // channel: [optional]  [ germline_resource ] or []
    germline_resource_tbi      // channel: [optional]  [ germline_resource_tbi ] or []
    germline_resource_provided // boolean: [mandatory] whether a germline resource was provided, so contamination tables are expected
    joint                      // boolean: [mandatory] joint calling, calls are per patient

    main:
    LEARNREADORIENTATIONMODEL(calls.map { meta, _vcf, _tbi, _stats, f1r2 -> [meta, f1r2] })

    vcf_with_priors = calls
        .map { meta, vcf, tbi, stats, _f1r2 -> [meta, vcf, tbi, stats] }
        .join(LEARNREADORIENTATIONMODEL.out.artifactprior, failOnDuplicate: true, failOnMismatch: true)

    pileup_table_normal = channel.empty()
    pileup_table_tumor  = channel.empty()
    contamination_table = channel.empty()
    segmentation_table  = channel.empty()

    if (germline_resource_provided) {
        cram_intervals = cram
            .combine(intervals)
            .map { meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals_, num_intervals ->
                [meta + [num_intervals: num_intervals], normal_cram, normal_crai, tumor_cram, tumor_crai, intervals_]
            }

        // In joint mode meta.id may be the patient, so the pileup ids come from the sample keys.
        // Dropping the tumour-specific keys (tumor_id, contamination) lets a normal shared by several tumours be summarised once.
        pileup = cram_intervals.multiMap { meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals_ ->
            tumor: [meta + [id: meta.tumor_id ?: meta.sample], tumor_cram, tumor_crai, intervals_]
            normal: [meta - meta.subMap('tumor_id', 'contamination') + [id: meta.normal_id], normal_cram, normal_crai, intervals_]
        }

        pileup_normal = pileup.normal.filter { meta, _cram, _crai, _intervals -> meta.normal_id }.unique()

        GETPILEUPSUMMARIES_NORMAL(pileup_normal, fasta, fai, dict, germline_resource, germline_resource_tbi)
        GETPILEUPSUMMARIES_TUMOR(pileup.tumor, fasta, fai, dict, germline_resource, germline_resource_tbi)

        // Figuring out if there is one or more table(s) from the same sample
        pileup_table_normal_branch = GETPILEUPSUMMARIES_NORMAL.out.table.branch { meta, _table ->
            intervals: meta.num_intervals > 1
            no_intervals: meta.num_intervals <= 1
        }
        pileup_table_tumor_branch = GETPILEUPSUMMARIES_TUMOR.out.table.branch { meta, _table ->
            intervals: meta.num_intervals > 1
            no_intervals: meta.num_intervals <= 1
        }

        // Only when using intervals
        pileup_table_normal_to_merge = pileup_table_normal_branch.intervals.map { meta, table -> [groupKey(meta, meta.num_intervals), table] }.groupTuple()
        pileup_table_tumor_to_merge = pileup_table_tumor_branch.intervals.map { meta, table -> [groupKey(meta, meta.num_intervals), table] }.groupTuple()

        GATHERPILEUPSUMMARIES_NORMAL(pileup_table_normal_to_merge, dict.map { _meta, dict_ -> [dict_] })
        GATHERPILEUPSUMMARIES_TUMOR(pileup_table_tumor_to_merge, dict.map { _meta, dict_ -> [dict_] })

        pileup_table_normal = channel.empty()
            .mix(GATHERPILEUPSUMMARIES_NORMAL.out.table, pileup_table_normal_branch.no_intervals)
            .map { meta, table -> [meta - meta.subMap('num_intervals'), table] }
        pileup_table_tumor = channel.empty()
            .mix(GATHERPILEUPSUMMARIES_TUMOR.out.table, pileup_table_tumor_branch.no_intervals)
            .map { meta, table -> [meta - meta.subMap('num_intervals'), table] }

        pileup_table_tumor_branch_status = pileup_table_tumor.branch { meta, _table ->
            paired: meta.normal_id
            tumor_only: true
        }

        // One normal table per [patient, normal_id] is combined with every tumour of that normal.
        // The pair meta is the cram meta with the pair id, so it matches the calls meta outside joint mode.
        pileup_table_normal_by_normal = pileup_table_normal.map { meta, table -> [[meta.patient, meta.normal_id], table] }

        calculatecontamination_in = pileup_table_tumor_branch_status.paired
            .map { meta, table -> [[meta.patient, meta.normal_id], meta + [id: "${meta.tumor_id}_vs_${meta.normal_id}".toString()], table] }
            .combine(pileup_table_normal_by_normal, by: 0)
            .map { _key, meta, tumor_table, normal_table -> [meta, tumor_table, normal_table] }
            .mix(pileup_table_tumor_branch_status.tumor_only.map { meta, table -> [meta, table, []] })

        CALCULATECONTAMINATION(calculatecontamination_in)

        if (joint) {
            // One table per input row, so each patient's group is released as soon as its own tables are done
            tables_per_patient = cram
                .map { meta, _normal_cram, _normal_crai, _tumor_cram, _tumor_crai -> [meta.patient, meta] }
                .groupTuple()
                .map { patient, metas -> [patient, metas.size()] }

            contamination_table = groupTablesByPatient(CALCULATECONTAMINATION.out.contamination, tables_per_patient)
            segmentation_table = groupTablesByPatient(CALCULATECONTAMINATION.out.segmentation, tables_per_patient)
        }
        else {
            contamination_table = CALCULATECONTAMINATION.out.contamination
            segmentation_table = CALCULATECONTAMINATION.out.segmentation
        }

        vcf_to_filter = vcf_with_priors
            .join(segmentation_table, failOnDuplicate: true, failOnMismatch: true)
            .join(contamination_table, failOnDuplicate: true, failOnMismatch: true)
            .map { meta, vcf, tbi, stats, artifactprior, segmentation, contamination ->
                [meta, vcf, tbi, stats, artifactprior, segmentation, contamination, []]
            }
    }
    else {
        // No contamination tables to wait for, so no join
        vcf_to_filter = vcf_with_priors.map { meta, vcf, tbi, stats, artifactprior ->
            [meta, vcf, tbi, stats, artifactprior, [], [], []]
        }
    }

    FILTERMUTECTCALLS(vcf_to_filter, fasta, fai, dict)

    emit:
    vcf                 = FILTERMUTECTCALLS.out.vcf                     // channel: [ meta, vcf ]
    tbi                 = FILTERMUTECTCALLS.out.tbi                     // channel: [ meta, tbi ]
    stats_filtered      = FILTERMUTECTCALLS.out.stats                   // channel: [ meta, stats ]
    artifact_priors     = LEARNREADORIENTATIONMODEL.out.artifactprior   // channel: [ meta, artifactprior ]
    contamination_table                                                 // channel: [ meta, contamination ]
    segmentation_table                                                  // channel: [ meta, segmentation ]
    pileup_table_normal                                                 // channel: [ meta, table ]
    pileup_table_tumor                                                  // channel: [ meta, table ]
}

// Joint calls are per patient: drop the per-sample keys as the caller does and group the tables by patient
def groupTablesByPatient(tables, tables_per_patient) {
    return tables
        .map { meta, table -> [meta.patient, meta - meta.subMap('normal_id', 'sample', 'tumor_id') + [id: meta.patient], table] }
        .combine(tables_per_patient, by: 0)
        .map { _patient, meta, table, n -> [groupKey(meta, n), table] }
        .groupTuple()
        .map { key, grouped_tables -> [key.getGroupTarget(), grouped_tables] }
}
