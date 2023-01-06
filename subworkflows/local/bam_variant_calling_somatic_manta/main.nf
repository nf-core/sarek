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

    versions = Channel.empty()

    MANTA_SOMATIC(cram, fasta, fasta_fai)

    // Figure out if using intervals or no_intervals
    candidate_small_indels_vcf = MANTA_SOMATIC.out.candidate_small_indels_vcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    candidate_small_indels_vcf_tbi = MANTA_SOMATIC.out.candidate_small_indels_vcf_tbi.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    manta_candidate_sv_vcf = MANTA_SOMATIC.out.candidate_sv_vcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    manta_diploid_sv_vcf = MANTA_SOMATIC.out.diploid_sv_vcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    manta_somatic_sv_vcf = MANTA_SOMATIC.out.somatic_sv_vcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    //Only when using intervals
    MERGE_MANTA_SV(
        candidate_small_indels_vcf.intervals.map{ meta, vcf ->
            [ groupKey(
                meta.subMap('normal_id', 'num_intervals', 'patient', 'sample', 'sex', 'tumor_id')
                    + [ id:meta.tumor_id + "_vs_" + meta.normal_id ],
                        meta.num_intervals),
            vcf ]
        }.groupTuple(),
        dict.map{ it -> [ [ id:'dict' ], it ] })

    MERGE_MANTA_SMALL_INDELS(
        manta_candidate_sv_vcf.intervals.map{ meta, vcf ->
            [ groupKey(
                meta.subMap('normal_id', 'num_intervals', 'patient', 'sample', 'sex', 'tumor_id')
                    + [ id:meta.tumor_id + "_vs_" + meta.normal_id ],
                        meta.num_intervals),
            vcf ]
        }.groupTuple(),
        dict.map{ it -> [ [ id:'dict' ], it ] })

    MERGE_MANTA_DIPLOID(
        manta_diploid_sv_vcf.intervals.map{ meta, vcf ->
            [ groupKey(
                meta.subMap('normal_id', 'num_intervals', 'patient', 'sample', 'sex', 'tumor_id')
                    + [ id:meta.tumor_id + "_vs_" + meta.normal_id ],
                        meta.num_intervals),
            vcf ]
        }.groupTuple(),
        dict.map{ it -> [ [ id:'dict' ], it ] })

    MERGE_MANTA_SOMATIC(
        manta_somatic_sv_vcf.intervals.map{ meta, vcf ->
            [ groupKey(
                meta.subMap('normal_id', 'num_intervals', 'patient', 'sample', 'sex', 'tumor_id')
                    + [ id:meta.tumor_id + "_vs_" + meta.normal_id ],
                        meta.num_intervals),
            vcf ]
        }.groupTuple(),
        dict.map{ it -> [ [ id:'dict' ], it ] })

    // Mix output channels for "no intervals" and "with intervals" results

    vcf = Channel.empty().mix(
        MERGE_MANTA_DIPLOID.out.vcf,
        MERGE_MANTA_SOMATIC.out.vcf,
        manta_diploid_sv_vcf.no_intervals,
        manta_somatic_sv_vcf.no_intervals
    ).map{ meta, vcf ->
            [ meta.subMap('normal_id', 'num_intervals', 'patient', 'sex', 'tumor_id')
                + [ id:meta.tumor_id + "_vs_" + meta.normal_id, variantcaller: 'manta' ],
                vcf ]
    }

    // Don't set variantcaller & num_intervals key. These files are not annotated, so they don't need it and joining with reads for StrelkaBP then fails
    candidate_small_indels_vcf = Channel.empty().mix(
        MERGE_MANTA_SMALL_INDELS.out.vcf,
        candidate_small_indels_vcf.no_intervals
    ).map{ meta, vcf ->
            [ meta.subMap('normal_id', 'patient', 'sex', 'tumor_id')
                + [ id:meta.tumor_id + "_vs_" + meta.normal_id ],
                vcf ]
    }

    candidate_small_indels_vcf_tbi = Channel.empty().mix(
        MERGE_MANTA_SMALL_INDELS.out.tbi,
        candidate_small_indels_vcf_tbi.no_intervals
    ).map{ meta, vcf ->
            [ meta.subMap('normal_id', 'patient', 'sex', 'tumor_id')
                + [ id:meta.tumor_id + "_vs_" + meta.normal_id ],
                vcf ]
    }

    versions = versions.mix(MERGE_MANTA_SV.out.versions)
    versions = versions.mix(MERGE_MANTA_SMALL_INDELS.out.versions)
    versions = versions.mix(MERGE_MANTA_DIPLOID.out.versions)
    versions = versions.mix(MERGE_MANTA_SOMATIC.out.versions)
    versions = versions.mix(MANTA_SOMATIC.out.versions)

    emit:
    candidate_small_indels_vcf
    candidate_small_indels_vcf_tbi
    vcf

    versions
}
