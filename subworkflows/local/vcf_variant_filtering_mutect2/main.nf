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

workflow VCF_VARIANT_FILTERING_MUTECT2 {
    take:
    calls_in                        // channel: [mandatory] [ meta, vcf, tbi, stats, f1r2 ] one row per pair, tumour-only sample or patient (joint)
    cram_in                         // channel: [mandatory] [ meta, normal_cram, normal_crai, tumour_cram, tumour_crai ] one row per pair or tumour-only sample ([] for the normal)
    intervals                       // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals
    fasta                           // channel: [mandatory] [ meta, fasta ]
    fai                             // channel: [mandatory] [ meta, fai ]
    dict                            // channel: [mandatory] [ meta, dict ]
    contamination_resource          // channel: [optional]  [ contamination_resource ] or [], common biallelic SNPs for getpileupsummaries
    contamination_resource_tbi      // channel: [optional]  [ contamination_resource_tbi ] or []
    contamination_resource_provided // boolean: [mandatory] whether a contamination resource was provided, so contamination tables are expected
    joint                           // boolean: [mandatory] joint calling, calls are per patient
    variantcaller                   // string:  [mandatory] caller of the calls, set as meta.variantcaller on all outputs

    main:
    calls = calls_in.map { meta, vcf, tbi, stats, f1r2 -> [meta + [variantcaller: variantcaller], vcf, tbi, stats, f1r2] }
    cram = cram_in.map { meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> [meta + [variantcaller: variantcaller], normal_cram, normal_crai, tumor_cram, tumor_crai] }

    LEARNREADORIENTATIONMODEL(calls.map { meta, _vcf, _tbi, _stats, f1r2 -> [meta, f1r2] })

    vcf_with_priors = calls
        .map { meta, vcf, tbi, stats, _f1r2 -> [meta, vcf, tbi, stats] }
        .join(LEARNREADORIENTATIONMODEL.out.artifactprior, failOnDuplicate: true, failOnMismatch: true)

    pileup_table_normal = channel.empty()
    pileup_table_tumor  = channel.empty()
    contamination_table = channel.empty()
    segmentation_table  = channel.empty()

    if (contamination_resource_provided) {
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

        // Tumour-only rows have no normal to summarise
        pileup_normal = pileup.normal.filter { meta, _cram, _crai, _intervals -> meta.normal_id }.unique()

        GETPILEUPSUMMARIES_NORMAL(pileup_normal, fasta, fai, dict, contamination_resource, contamination_resource_tbi)
        GETPILEUPSUMMARIES_TUMOR(pileup.tumor, fasta, fai, dict, contamination_resource, contamination_resource_tbi)

        pileup_normal_gather = gatherPileups(GETPILEUPSUMMARIES_NORMAL.out.table)
        pileup_tumor_gather = gatherPileups(GETPILEUPSUMMARIES_TUMOR.out.table)

        GATHERPILEUPSUMMARIES_NORMAL(pileup_normal_gather.intervals, dict.map { _meta, dict_ -> [dict_] })
        GATHERPILEUPSUMMARIES_TUMOR(pileup_tumor_gather.intervals, dict.map { _meta, dict_ -> [dict_] })

        pileup_table_normal = channel.empty()
            .mix(GATHERPILEUPSUMMARIES_NORMAL.out.table, pileup_normal_gather.no_intervals)
            .map { meta, table -> [meta - meta.subMap('num_intervals'), table] }
        pileup_table_tumor = channel.empty()
            .mix(GATHERPILEUPSUMMARIES_TUMOR.out.table, pileup_tumor_gather.no_intervals)
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
            // Joint calls are per patient, so are their tables. Without a group size, groupTuple waits for the tables of all patients.
            contamination_table = CALCULATECONTAMINATION.out.contamination.map { meta, table -> [patientMeta(meta), table] }.groupTuple()
            segmentation_table = CALCULATECONTAMINATION.out.segmentation.map { meta, table -> [patientMeta(meta), table] }.groupTuple()
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

    fasta_fai_dict = fasta
        .combine(fai)
        .combine(dict)
        .map { meta, fasta_, _meta_fai, fai_, _meta_dict, dict_ -> [meta, fasta_, fai_, dict_] }
        .collect()

    FILTERMUTECTCALLS(vcf_to_filter, fasta_fai_dict)

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

// Patient-level meta of joint calling: the per-sample keys are dropped and the id is the patient.
// Used for the joint Mutect2 calls and for their contamination tables, which are joined on it.
// data_type (bam or cram) is per sample in tumour-only meta: samples of one patient may differ.
def patientMeta(meta) {
    return meta - meta.subMap('normal_id', 'tumor_id', 'sample', 'contamination', 'data_type') + [id: meta.patient]
}

// One row per call or sample with the outputs of all its intervals (a single one without intervals, num_intervals 0 or 1).
// groupKey releases each row as soon as its intervals are done.
def groupByIntervals(ch) {
    return ch
        .map { meta, file -> [groupKey(meta, Math.max(meta.num_intervals, 1)), file] }
        .groupTuple()
}

// Scattered pileup tables of a sample are to be gathered, a single table passes through
def gatherPileups(tables) {
    return groupByIntervals(tables).branch { meta, tables_ ->
        intervals: meta.num_intervals > 1
        no_intervals: true
            return [meta, tables_[0]]
    }
}
