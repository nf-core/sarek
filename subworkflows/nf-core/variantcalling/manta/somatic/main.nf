include { GATK4_MERGEVCFS as MERGE_MANTA_DIPLOID           } from '../../../../../modules/nf-core/modules/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_MANTA_SMALL_INDELS      } from '../../../../../modules/nf-core/modules/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_MANTA_SOMATIC           } from '../../../../../modules/nf-core/modules/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_MANTA_SV                } from '../../../../../modules/nf-core/modules/gatk4/mergevcfs/main'
include { MANTA_SOMATIC                                    } from '../../../../../modules/nf-core/modules/manta/somatic/main'

workflow RUN_MANTA_SOMATIC {
    take:
    cram                     // channel: [mandatory] [meta, normal_cram, normal_crai, tumor_cram, tumor_crai, interval.bed.gz, interval.bed.gz.tbi]
    dict                     // channel: [optional]
    fasta                    // channel: [mandatory]
    fasta_fai                // channel: [mandatory]
    intervals_bed_gz         // channel: [optional]  Contains a bed.gz file of all intervals combined provided with the cram input(s). Mandatory if interval files are used.

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
<<<<<<< HEAD
    BGZIP_VC_MANTA_SMALL_INDELS(manta_candidate_small_indels_vcf.intervals)

    CONCAT_MANTA_SMALL_INDELS(
        BGZIP_VC_MANTA_SMALL_INDELS.out.output.map{ meta, vcf ->
=======
    MERGE_MANTA_SV(
        manta_candidate_small_indels_vcf.intervals.map{ meta, vcf ->
>>>>>>> origin/dev

                [groupKey([patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, gender:meta.gender, id:meta.tumor_id + "_vs_" + meta.normal_id, num_intervals:meta.num_intervals],
                            meta.num_intervals),
                vcf]

            }.groupTuple(),
        dict)

<<<<<<< HEAD
    BGZIP_VC_MANTA_SV(manta_candidate_sv_vcf.intervals)

    CONCAT_MANTA_SV(
        BGZIP_VC_MANTA_SV.out.output.map{ meta, vcf ->

                [groupKey([patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, gender:meta.gender, id:meta.tumor_id + "_vs_" + meta.normal_id, num_intervals:meta.num_intervals],
                            meta.num_intervals),
                vcf]
=======
    MERGE_MANTA_SMALL_INDELS(
        manta_candidate_sv_vcf.intervals.map{ meta, vcf ->
                new_meta = [patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, gender:meta.gender, id:meta.tumor_id + "_vs_" + meta.normal_id, num_intervals:meta.num_intervals]
>>>>>>> origin/dev

            }.groupTuple(),
        dict)

<<<<<<< HEAD
    BGZIP_VC_MANTA_DIPLOID(manta_diploid_sv_vcf.intervals)

    CONCAT_MANTA_DIPLOID(
        BGZIP_VC_MANTA_DIPLOID.out.output.map{ meta, vcf ->
=======
    MERGE_MANTA_DIPLOID(
        manta_diploid_sv_vcf.intervals.map{ meta, vcf ->
                new_meta = [patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, gender:meta.gender, id:meta.tumor_id + "_vs_" + meta.normal_id, num_intervals:meta.num_intervals]
>>>>>>> origin/dev

                [groupKey([patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, gender:meta.gender, id:meta.tumor_id + "_vs_" + meta.normal_id, num_intervals:meta.num_intervals],
                            meta.num_intervals),
                vcf]

            }.groupTuple(),
        dict)

<<<<<<< HEAD
    BGZIP_VC_MANTA_SOMATIC(manta_somatic_sv_vcf.intervals)

    CONCAT_MANTA_SOMATIC(
        BGZIP_VC_MANTA_SOMATIC.out.output.map{ meta, vcf ->
=======
    MERGE_MANTA_SOMATIC(
        manta_somatic_sv_vcf.intervals.map{ meta, vcf ->
                new_meta = [patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, gender:meta.gender, id:meta.tumor_id + "_vs_" + meta.normal_id, num_intervals:meta.num_intervals]
>>>>>>> origin/dev

                [groupKey([patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, gender:meta.gender, id:meta.tumor_id + "_vs_" + meta.normal_id, num_intervals:meta.num_intervals],
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
        [[patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, gender:meta.gender, id:meta.tumor_id + "_vs_" + meta.normal_id, num_intervals:meta.num_intervals, variantcaller:"Manta"],
        vcf]
    }

    // Don't set variantcaller & num_intervals key. These files are not annotated, so they don't need it and joining with reads for StrelkaBP then fails
    manta_candidate_small_indels_vcf = Channel.empty().mix(
        MERGE_MANTA_SMALL_INDELS.out.vcf,
        manta_candidate_small_indels_vcf.no_intervals
    ).map{ meta, vcf ->
        [[patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, gender:meta.gender, id:meta.tumor_id + "_vs_" + meta.normal_id],
        vcf]
    }

    manta_candidate_small_indels_vcf_tbi = Channel.empty().mix(
        MERGE_MANTA_SMALL_INDELS.out.tbi,
        manta_candidate_small_indels_vcf_tbi.no_intervals
    ).map{ meta, vcf ->
        [[patient:meta.patient, normal_id:meta.normal_id, tumor_id:meta.tumor_id, gender:meta.gender, id:meta.tumor_id + "_vs_" + meta.normal_id],
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
