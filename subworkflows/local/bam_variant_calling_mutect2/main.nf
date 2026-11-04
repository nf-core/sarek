//
// MUTECT2: tumour-normal and tumour-only variant calling, filtering in VCF_VARIANT_FILTERING_MUTECT2
//

include { VCF_VARIANT_FILTERING_MUTECT2              } from '../vcf_variant_filtering_mutect2'
include { GATK4_MERGEMUTECTSTATS as MERGEMUTECTSTATS } from '../../../modules/nf-core/gatk4/mergemutectstats'
include { GATK4_MERGEVCFS as MERGE_MUTECT2           } from '../../../modules/nf-core/gatk4/mergevcfs'
include { GATK4_MUTECT2 as MUTECT2                   } from '../../../modules/nf-core/gatk4/mutect2'

workflow BAM_VARIANT_CALLING_MUTECT2 {
    take:
    cram                            // channel: [mandatory] [ meta, normal_cram, normal_crai, tumor_cram, tumor_crai ] one row per pair or tumour-only sample ([] for the normal)
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
        // One call per patient; a normal shared by several tumours is passed once.
        // normal_id is kept for --normal-sample and dropped after calling.
        input = cram
            .flatMap { meta, normal_cram, normal_crai, tumor_cram, tumor_crai ->
                def patient_meta = meta - meta.subMap('normal_id', 'tumor_id', 'sample', 'contamination') + [id: meta.patient] + (meta.normal_id ? [normal_id: meta.normal_id] : [:])
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

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_branch = MUTECT2.out.vcf.branch { meta, _vcf ->
        intervals: meta.num_intervals > 1
        no_intervals: meta.num_intervals <= 1
    }

    // Figuring out if there is one or more tbi(s) from the same sample
    tbi_branch = MUTECT2.out.tbi.branch { meta, _tbi ->
        intervals: meta.num_intervals > 1
        no_intervals: meta.num_intervals <= 1
    }

    // Figuring out if there is one or more stats(s) from the same sample
    stats_branch = MUTECT2.out.stats.branch { meta, _stats ->
        intervals: meta.num_intervals > 1
        no_intervals: meta.num_intervals <= 1
    }

    // Figuring out if there is one or more f1r2(s) from the same sample
    f1r2_branch = MUTECT2.out.f1r2.branch { meta, _f1r2 ->
        intervals: meta.num_intervals > 1
        no_intervals: meta.num_intervals <= 1
    }

    // Only when using intervals
    vcf_to_merge = vcf_branch.intervals.map { meta, vcf -> [groupKey(meta, meta.num_intervals), vcf] }.groupTuple()
    stats_to_merge = stats_branch.intervals.map { meta, stats -> [groupKey(meta, meta.num_intervals), stats] }.groupTuple()
    f1r2_to_merge = f1r2_branch.intervals.map { meta, f1r2 -> [groupKey(meta, meta.num_intervals), f1r2] }.groupTuple()

    MERGE_MUTECT2(vcf_to_merge, dict)
    MERGEMUTECTSTATS(stats_to_merge)

    // Mutect2 calls filtered by filtermutectcalls in the shared filtering subworkflow
    // meta joint calling:  [id:patient, patient, sex, variantcaller]
    // meta paired calling: [id:tumor_id_vs_normal_id, normal_id, patient, sex, tumor_id, variantcaller]
    // meta tumour-only:    [id:sample, patient, sample, sex, status, variantcaller]
    calls = channel.empty()
        .mix(MERGE_MUTECT2.out.vcf, vcf_branch.no_intervals)
        .join(channel.empty().mix(MERGE_MUTECT2.out.tbi, tbi_branch.no_intervals), failOnDuplicate: true, failOnMismatch: true)
        .join(channel.empty().mix(MERGEMUTECTSTATS.out.stats, stats_branch.no_intervals), failOnDuplicate: true, failOnMismatch: true)
        .join(channel.empty().mix(f1r2_to_merge, f1r2_branch.no_intervals), failOnDuplicate: true, failOnMismatch: true)
        .map { meta, vcf, tbi, stats, f1r2 ->
            [meta - meta.subMap(joint_mutect2 ? ['normal_id', 'num_intervals'] : ['num_intervals']) + [variantcaller: 'mutect2'], vcf, tbi, stats, f1r2]
        }

    VCF_VARIANT_FILTERING_MUTECT2(
        calls,
        cram.map { meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> [meta + [variantcaller: 'mutect2'], normal_cram, normal_crai, tumor_cram, tumor_crai] },
        intervals,
        fasta,
        fai,
        dict,
        contamination_resource,
        contamination_resource_tbi,
        contamination_resource_provided,
        joint_mutect2,
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
