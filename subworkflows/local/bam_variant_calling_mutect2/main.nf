//
// MUTECT2: tumor-normal and tumor-only variant calling, filtering in VCF_VARIANT_FILTERING_MUTECT2
//

include { VCF_VARIANT_FILTERING_MUTECT2              } from '../vcf_variant_filtering_mutect2'
include { groupByIntervals                           } from '../vcf_variant_filtering_mutect2'
include { patientMeta                                } from '../vcf_variant_filtering_mutect2'
include { GATK4_MERGEMUTECTSTATS as MERGEMUTECTSTATS } from '../../../modules/nf-core/gatk4/mergemutectstats'
include { GATK4_MERGEVCFS as MERGE_MUTECT2           } from '../../../modules/nf-core/gatk4/mergevcfs'
include { GATK4_MUTECT2 as MUTECT2                   } from '../../../modules/nf-core/gatk4/mutect2'

workflow BAM_VARIANT_CALLING_MUTECT2 {
    take:
    cram                            // channel: [mandatory] [ meta, normal_cram, normal_crai, tumor_cram, tumor_crai ] one row per pair or tumor-only sample ([] for the normal)
    fasta                           // channel: [mandatory] [ meta, fasta ]
    fai                             // channel: [mandatory] [ meta, fai ]
    dict                            // channel: [mandatory] [ meta, dict ]
    force_call_alleles              // channel: [optional]  [ force_call_alleles ] or []
    force_call_alleles_tbi          // channel: [optional]  [ force_call_alleles_tbi ] or []
    germline_resource               // channel: [optional]  [ germline_resource ] or []
    germline_resource_tbi           // channel: [optional]  [ germline_resource_tbi ] or []
    contamination_resource          // channel: [optional]  [ contamination_resource ] or [], common biallelic SNPs for getpileupsummaries
    contamination_resource_tbi      // channel: [optional]  [ contamination_resource_tbi ] or []
    panel_of_normals                // channel: [optional]  [ panel_of_normals ] or []
    panel_of_normals_tbi            // channel: [optional]  [ panel_of_normals_tbi ] or []
    intervals                       // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals
    joint_mutect2                   // boolean: [mandatory] call all samples of a patient together
    contamination_resource_provided // boolean: [mandatory] whether a contamination resource was provided, so contamination tables are expected

    main:
    if (joint_mutect2) {
        // One call per patient; a normal shared by several tumors is passed once.
        // normal_id is kept for --normal-sample and dropped after calling.
        input = cram
            .flatMap { meta, normal_cram, normal_crai, tumor_cram, tumor_crai ->
                def patient_meta = patientMeta(meta) + (meta.normal_id ? [normal_id: meta.normal_id] : [:])
                normal_cram ? [[patient_meta, normal_cram, normal_crai], [patient_meta, tumor_cram, tumor_crai]] : [[patient_meta, tumor_cram, tumor_crai]]
            }
            .unique()
            .groupTuple()
    }
    else {
        input = cram.map { meta, normal_cram, normal_crai, tumor_cram, tumor_crai ->
            [meta, [normal_cram, tumor_cram].findAll(), [normal_crai, tumor_crai].findAll()]
        }
    }

    // Move num_intervals to meta map for scatter and gather
    input_intervals = input
        .combine(intervals)
        .map { meta, input_list, input_index_list, intervals_, num_intervals -> [meta + [num_intervals: num_intervals], input_list, input_index_list, intervals_] }

    MUTECT2(input_intervals, fasta, fai.map { meta, index -> [ meta, index, [] ] }, dict, force_call_alleles, force_call_alleles_tbi, germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi)

    // All outputs of a call, scattered or not, in one row
    mutect2_out = groupByIntervals(MUTECT2.out.vcf)
        .join(groupByIntervals(MUTECT2.out.tbi), failOnDuplicate: true, failOnMismatch: true)
        .join(groupByIntervals(MUTECT2.out.stats), failOnDuplicate: true, failOnMismatch: true)
        .join(groupByIntervals(MUTECT2.out.f1r2), failOnDuplicate: true, failOnMismatch: true)
        .branch { meta, vcf, tbi, stats, f1r2 ->
            // Scattered: vcf and stats are merged (MergeVcfs writes the tbi), all f1r2 go to LearnReadOrientationModel
            intervals: meta.num_intervals > 1
                return [meta, vcf, stats, f1r2]
            no_intervals: true
                return [meta, vcf[0], tbi[0], stats[0], f1r2[0]]
        }

    MERGE_MUTECT2(mutect2_out.intervals.map { meta, vcf, _stats, _f1r2 -> [meta, vcf] }, dict)
    MERGEMUTECTSTATS(mutect2_out.intervals.map { meta, _vcf, stats, _f1r2 -> [meta, stats] })

    // Mutect2 calls filtered by filtermutectcalls in the shared filtering subworkflow, which sets meta.variantcaller
    // meta joint calling:  [id:patient, patient, sex]
    // meta paired calling: [id:tumor_id_vs_normal_id, normal_id, patient, sex, tumor_id]
    // meta tumor-only:    [id:sample, data_type, patient, sample, sex, status]
    calls = MERGE_MUTECT2.out.vcf
        .join(MERGE_MUTECT2.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .join(MERGEMUTECTSTATS.out.stats, failOnDuplicate: true, failOnMismatch: true)
        .join(mutect2_out.intervals.map { meta, _vcf, _stats, f1r2 -> [meta, f1r2] }, failOnDuplicate: true, failOnMismatch: true)
        .mix(mutect2_out.no_intervals)
        .map { meta, vcf, tbi, stats, f1r2 ->
            [meta - meta.subMap(joint_mutect2 ? ['normal_id', 'num_intervals'] : ['num_intervals']), vcf, tbi, stats, f1r2]
        }

    VCF_VARIANT_FILTERING_MUTECT2(
        calls,
        cram,
        intervals,
        fasta,
        fai,
        dict,
        contamination_resource,
        contamination_resource_tbi,
        contamination_resource_provided,
        joint_mutect2,
        'mutect2',
    )

    emit:
    vcf                 = VCF_VARIANT_FILTERING_MUTECT2.out.vcf                 // channel: [ meta, vcf ]
    tbi                 = VCF_VARIANT_FILTERING_MUTECT2.out.tbi                 // channel: [ meta, tbi ]
    stats_filtered      = VCF_VARIANT_FILTERING_MUTECT2.out.stats_filtered      // channel: [ meta, stats ]
    artifact_priors     = VCF_VARIANT_FILTERING_MUTECT2.out.artifact_priors     // channel: [ meta, artifactprior ]
    pileup_table_normal = VCF_VARIANT_FILTERING_MUTECT2.out.pileup_table_normal // channel: [ meta, table ]
    pileup_table_tumor  = VCF_VARIANT_FILTERING_MUTECT2.out.pileup_table_tumor  // channel: [ meta, table ]
    contamination_table = VCF_VARIANT_FILTERING_MUTECT2.out.contamination_table // channel: [ meta, contamination ]
    segmentation_table  = VCF_VARIANT_FILTERING_MUTECT2.out.segmentation_table  // channel: [ meta, segmentation ]
}
