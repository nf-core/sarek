include { GATK4_MERGEVCFS as MERGE_MANTA_DIPLOID           } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_MANTA_SMALL_INDELS      } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_MANTA_SOMATIC           } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_MANTA_SV                } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { MANTA_SOMATIC                                    } from '../../../modules/nf-core/manta/somatic/main'

workflow BAM_VARIANT_CALLING_SOMATIC_MANTA {
    take:
    cram                     // channel: [mandatory] [meta, normal_cram, normal_crai, tumor_cram, tumor_crai, interval.bed.gz, interval.bed.gz.tbi]
    dict                     // channel: [optional]
    fasta                    // channel: [mandatory]
    fasta_fai                // channel: [mandatory]

    main:

    ch_versions = Channel.empty()

    MANTA_SOMATIC(cram, fasta, fasta_fai)

    // Figure out if using intervals or no_intervals
    MANTA_SOMATIC.out.candidate_small_indels_vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{manta_candidate_small_indels_vcf}

    MANTA_SOMATIC.out.candidate_small_indels_vcf_tbi.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{manta_candidate_small_indels_vcf_tbi}

    MANTA_SOMATIC.out.candidate_sv_vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{manta_candidate_sv_vcf}

    MANTA_SOMATIC.out.diploid_sv_vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{manta_diploid_sv_vcf}

    MANTA_SOMATIC.out.somatic_sv_vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{manta_somatic_sv_vcf}

    //Only when using intervals
    MERGE_MANTA_SV(
        manta_candidate_small_indels_vcf.intervals.map{ meta, vcf ->

                [groupKey([
                            id:             meta.tumor_id + "_vs_" + meta.normal_id,
                            normal_id:      meta.normal_id,
                            num_intervals:  meta.num_intervals,
                            patient:        meta.patient,
                            sex:            meta.sex,
                            tumor_id:       meta.tumor_id
                        ],
                        meta.num_intervals),
                vcf]

            }.groupTuple(),
        dict)

    MERGE_MANTA_SMALL_INDELS(
        manta_candidate_sv_vcf.intervals.map{ meta, vcf ->

                [groupKey([
                            id:             meta.tumor_id + "_vs_" + meta.normal_id,
                            normal_id:      meta.normal_id,
                            num_intervals:  meta.num_intervals,
                            patient:        meta.patient,
                            sex:            meta.sex,
                            tumor_id:       meta.tumor_id
                        ],
                        meta.num_intervals),
                vcf]

            }.groupTuple(),
        dict)

    MERGE_MANTA_DIPLOID(
        manta_diploid_sv_vcf.intervals.map{ meta, vcf ->
                new_meta = [
                            id:             meta.tumor_id + "_vs_" + meta.normal_id,
                            normal_id:      meta.normal_id,
                            num_intervals:  meta.num_intervals,
                            patient:        meta.patient,
                            sex:            meta.sex,
                            tumor_id:       meta.tumor_id
                        ]

                [groupKey([
                            id:             meta.tumor_id + "_vs_" + meta.normal_id,
                            normal_id:      meta.normal_id,
                            num_intervals:  meta.num_intervals,
                            patient:        meta.patient,
                            sex:            meta.sex,
                            tumor_id:       meta.tumor_id
                        ],
                        meta.num_intervals),
                vcf]

            }.groupTuple(),
        dict)

    MERGE_MANTA_SOMATIC(
        manta_somatic_sv_vcf.intervals.map{ meta, vcf ->

                [groupKey([
                            id:             meta.tumor_id + "_vs_" + meta.normal_id,
                            normal_id:      meta.normal_id,
                            num_intervals:  meta.num_intervals,
                            patient:        meta.patient,
                            sex:            meta.sex,
                            tumor_id:       meta.tumor_id
                        ],
                        meta.num_intervals),
                vcf]

            }.groupTuple(),
        dict)

    // Mix output channels for "no intervals" and "with intervals" results
    manta_vcf = Channel.empty().mix(
        MERGE_MANTA_DIPLOID.out.vcf,
        MERGE_MANTA_SOMATIC.out.vcf,
        manta_diploid_sv_vcf.no_intervals,
        manta_somatic_sv_vcf.no_intervals
    ).map{ meta, vcf ->
        [[
            id:             meta.tumor_id + "_vs_" + meta.normal_id,
            num_intervals:  meta.num_intervals,
            normal_id:      meta.normal_id,
            patient:        meta.patient,
            sex:            meta.sex,
            tumor_id:       meta.tumor_id,
            variantcaller:  "manta"
        ],
        vcf]
    }

    // Don't set variantcaller & num_intervals key. These files are not annotated, so they don't need it and joining with reads for StrelkaBP then fails
    manta_candidate_small_indels_vcf = Channel.empty().mix(
        MERGE_MANTA_SMALL_INDELS.out.vcf,
        manta_candidate_small_indels_vcf.no_intervals
    ).map{ meta, vcf ->
        [[
            id:         meta.tumor_id + "_vs_" + meta.normal_id,
            normal_id:  meta.normal_id,
            patient:    meta.patient,
            sex:        meta.sex,
            tumor_id:   meta.tumor_id,
        ],
        vcf]
    }

    manta_candidate_small_indels_vcf_tbi = Channel.empty().mix(
        MERGE_MANTA_SMALL_INDELS.out.tbi,
        manta_candidate_small_indels_vcf_tbi.no_intervals
    ).map{ meta, vcf ->
        [[
            id:         meta.tumor_id + "_vs_" + meta.normal_id,
            normal_id:  meta.normal_id,
            patient:    meta.patient,
            sex:        meta.sex,
            tumor_id:   meta.tumor_id
        ],
        vcf]
    }

    ch_versions = ch_versions.mix(MERGE_MANTA_SV.out.versions)
    ch_versions = ch_versions.mix(MERGE_MANTA_SMALL_INDELS.out.versions)
    ch_versions = ch_versions.mix(MERGE_MANTA_DIPLOID.out.versions)
    ch_versions = ch_versions.mix(MERGE_MANTA_SOMATIC.out.versions)
    ch_versions = ch_versions.mix(MANTA_SOMATIC.out.versions)

    emit:
    manta_vcf
    manta_candidate_small_indels_vcf
    manta_candidate_small_indels_vcf_tbi
    versions = ch_versions

}
