include { GATK4_MERGEVCFS as MERGE_MANTA_DIPLOID      } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_MANTA_SMALL_INDELS } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_MANTA_SV           } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { MANTA_GERMLINE                              } from '../../../modules/nf-core/manta/germline/main'

// Seems to be the consensus on upstream modules implementation too
workflow BAM_VARIANT_CALLING_GERMLINE_MANTA {
    take:
    cram                     // channel: [mandatory] [meta, cram, crai, interval.bed.gz, interval.bed.gz.tbi]
    dict                     // channel: [optional]
    fasta                    // channel: [mandatory]
    fasta_fai                // channel: [mandatory]

    main:

    ch_versions = Channel.empty()

    MANTA_GERMLINE(cram, fasta, fasta_fai)

    // Figure out if using intervals or no_intervals
    MANTA_GERMLINE.out.candidate_small_indels_vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{manta_small_indels_vcf}

    MANTA_GERMLINE.out.candidate_sv_vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{manta_sv_vcf}

    MANTA_GERMLINE.out.diploid_sv_vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{manta_diploid_sv_vcf}

    // Only when using intervals
    MERGE_MANTA_SMALL_INDELS(
        manta_small_indels_vcf.intervals
            .map{ meta, vcf ->

                [groupKey([
                            id:             meta.sample,
                            num_intervals:  meta.num_intervals,
                            patient:        meta.patient,
                            sample:         meta.sample,
                            sex:            meta.sex,
                            status:         meta.status,
                        ],
                        meta.num_intervals),
                vcf]
            }.groupTuple(),
        dict)

    MERGE_MANTA_SV(
        manta_sv_vcf.intervals
            .map{ meta, vcf ->

                [groupKey([
                            id:             meta.sample,
                            num_intervals:  meta.num_intervals,
                            patient:        meta.patient,
                            sample:         meta.sample,
                            sex:            meta.sex,
                            status:         meta.status,
                        ],
                        meta.num_intervals),
                vcf]

            }.groupTuple(),
        dict)

    MERGE_MANTA_DIPLOID(
        manta_diploid_sv_vcf.intervals
            .map{ meta, vcf ->

                [groupKey([
                            id:             meta.sample,
                            num_intervals:  meta.num_intervals,
                            patient:        meta.patient,
                            sample:         meta.sample,
                            status:         meta.status,
                            sex:            meta.sex,
                        ],
                        meta.num_intervals),
                vcf]

            }.groupTuple(),
        dict)

    // Mix output channels for "no intervals" and "with intervals" results
    // Only diploid SV should get annotated
    manta_vcf = Channel.empty().mix(
                    MERGE_MANTA_DIPLOID.out.vcf,
                    manta_diploid_sv_vcf.no_intervals)
                .map{ meta, vcf ->
                    [[
                        id:             meta.sample,
                        num_intervals:  meta.num_intervals,
                        patient:        meta.patient,
                        sample:         meta.sample,
                        status:         meta.status,
                        sex:            meta.sex,
                        variantcaller:  "manta"],
                    vcf]
                }

    ch_versions = ch_versions.mix(MERGE_MANTA_DIPLOID.out.versions)
    ch_versions = ch_versions.mix(MERGE_MANTA_SMALL_INDELS.out.versions)
    ch_versions = ch_versions.mix(MERGE_MANTA_SV.out.versions)
    ch_versions = ch_versions.mix(MANTA_GERMLINE.out.versions)

    emit:
    manta_vcf
    versions = ch_versions
}
