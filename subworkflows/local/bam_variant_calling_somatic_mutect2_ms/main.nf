//
// Run GATK mutect2 in multi sample mode, getepileupsummaries, calculatecontamination, learnreadorientationmodel and filtermutectcalls
//

include { GATK4_MUTECT2                   as MUTECT2_MS               } from '../../../modules/nf-core/gatk4/mutect2/main'
include { GATK4_MERGEVCFS                 as MERGE_MUTECT2            } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEMUTECTSTATS          as MERGEMUTECTSTATS         } from '../../../modules/nf-core/gatk4/mergemutectstats/main'

workflow BAM_VARIANT_CALLING_SOMATIC_MUTECT2_MS {
    take:
    input                     // channel: [ val(meta), normal_cram, normal_crai, tumor_cram, tumor_crai, intervals ]
    fasta                     // channel: /path/to/reference/fasta
    fai                       // channel: /path/to/reference/fasta/index
    dict                      // channel: /path/to/reference/fasta/dictionary
    germline_resource         // channel: /path/to/germline/resource
    germline_resource_tbi     // channel: /path/to/germline/index
    panel_of_normals          // channel: /path/to/panel/of/normals
    panel_of_normals_tbi      // channel: /path/to/panel/of/normals/index

    main:
    ch_versions = Channel.empty()

    // Prepare multi sample Mutect2 variant calling channels
    // Separate normal cram files
    ch_normal_cram = input.map{meta, n_cram, n_crai, t_cram, t_crai, intervals -> [[id: meta.patient, patient: meta.patient, normal_id: meta.normal_id, sex: meta.sex, num_intervals: meta.num_intervals], n_cram, n_crai]}.unique().groupTuple()
    // If there are multiple normal cram files, use the first one only
    ch_normal_cram_first = ch_normal_cram.map{it -> [it[0], it[1][0], it[2][0]]}
    // Separate tumor cram files
    ch_tumor_cram = input.map{meta, n_cram, n_crai, t_cram, t_crai, intervals -> [[id: meta.patient, patient: meta.patient, normal_id: meta.normal_id, sex: meta.sex, num_intervals: meta.num_intervals], t_cram, t_crai]}.unique()
    // Merge normal and tumor samples by patient
    ch_tn_cram = ch_normal_cram_first.mix(ch_tumor_cram).groupTuple()
    // Add intervals back
    ch_pt_intervals = input.map{meta, n_cram, n_crai, t_cram, t_crai, intervals -> [[id: meta.patient, patient: meta.patient, normal_id: meta.normal_id, sex: meta.sex, num_intervals: meta.num_intervals], intervals]}.unique().groupTuple()
    ch_tn_intervals = ch_tn_cram.join(ch_pt_intervals).transpose(by : [3])
    ch_tn_intervals.view()
 
    MUTECT2_MS(
            ch_tn_intervals,
            fasta,
            fai,
            dict,
            germline_resource,
            germline_resource_tbi,
            panel_of_normals,
            panel_of_normals_tbi)


    // Figure out if using intervals or no_intervals
    MUTECT2_MS.out.vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{ mutect2_vcf_branch }

    MUTECT2_MS.out.tbi.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{ mutect2_tbi_branch }

    MUTECT2_MS.out.stats.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{ mutect2_stats_branch }

    MUTECT2_MS.out.f1r2.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{ mutect2_f1r2_branch }

    //Only when using intervals
    MERGE_MUTECT2(
        mutect2_vcf_branch.intervals
        .map{ meta, vcf ->

            new_meta = [id:meta.patient,
                        normal_id:meta.normal_id,
                        num_intervals:meta.num_intervals,
                        patient:meta.patient,
                        sex:meta.sex
                    ]

            [groupKey(new_meta, meta.num_intervals), vcf]
        }.groupTuple(),
        dict
    )

    mutect2_vcf = Channel.empty().mix(
        MERGE_MUTECT2.out.vcf,
        mutect2_vcf_branch.no_intervals)

    mutect2_tbi = Channel.empty().mix(
        MERGE_MUTECT2.out.tbi,
        mutect2_tbi_branch.no_intervals)
    
    //Merge Mutect2 Stats
    MERGEMUTECTSTATS(
        mutect2_stats_branch.intervals
        .map{ meta, stats ->

            new_meta = [id:meta.patient,
                        normal_id:      meta.normal_id,
                        num_intervals:  meta.num_intervals,
                        patient:        meta.patient,
                        sex:            meta.sex
                    ]

            [groupKey(new_meta, meta.num_intervals), stats]
        }.groupTuple())

    mutect2_stats = Channel.empty().mix(
        MERGEMUTECTSTATS.out.stats,
        mutect2_stats_branch.no_intervals)


    ch_versions = ch_versions.mix(MUTECT2_MS.out.versions)
    ch_versions = ch_versions.mix(MERGE_MUTECT2.out.versions)
    ch_versions = ch_versions.mix(MERGEMUTECTSTATS.out.versions)

    emit:
    mutect2_vcf            = mutect2_vcf                                    // channel: [ val(meta), [ vcf ] ]
    mutect2_stats          = mutect2_stats                                  // channel: [ val(meta), [ stats ] ]

    versions               = ch_versions                                    // channel: [ versions.yml ]
}
