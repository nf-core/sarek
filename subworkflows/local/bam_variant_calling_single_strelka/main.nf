include { GATK4_MERGEVCFS as MERGE_STRELKA        } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_STRELKA_GENOME } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { STRELKA_GERMLINE as STRELKA_SINGLE      } from '../../../modules/nf-core/strelka/germline/main'

workflow BAM_VARIANT_CALLING_SINGLE_STRELKA {
    take:
    cram                     // channel: [mandatory] [meta, cram, crai, interval.bed.gz, interval.bed.gz.tbi]
    dict                     // channel: [optional]
    fasta                    // channel: [mandatory]
    fasta_fai                // channel: [mandatory]

    main:

    ch_versions = Channel.empty()

    STRELKA_SINGLE(cram, fasta, fasta_fai)

    // Figure out if using intervals or no_intervals
    STRELKA_SINGLE.out.vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{strelka_vcf}

    STRELKA_SINGLE.out.vcf_tbi.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{strelka_vcf_tbi}

    STRELKA_SINGLE.out.genome_vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{strelka_genome_vcf}

    MERGE_STRELKA(
        strelka_vcf.intervals
            .map{ meta, vcf ->
                new_meta = [
                                id:             meta.sample,
                                num_intervals:  meta.num_intervals,
                                patient:        meta.patient,
                                sample:         meta.sample,
                                sex:            meta.sex,
                                status:         meta.status
                            ]

                [groupKey(new_meta, meta.num_intervals), vcf]
            }.groupTuple(),
        dict
    )

    MERGE_STRELKA_GENOME(
        strelka_genome_vcf.intervals
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
        dict
    )

    // Mix output channels for "no intervals" and "with intervals" results
    // Only strelka variant vcf should get annotated
    strelka_vcf = Channel.empty().mix(
                    MERGE_STRELKA.out.vcf,
                    strelka_vcf.no_intervals)
                .map{ meta, vcf ->
                    [[
                        id:             meta.sample,
                        num_intervals:  meta.num_intervals,
                        patient:        meta.patient,
                        sample:         meta.sample,
                        sex:            meta.sex,
                        status:         meta.status,
                        variantcaller:  "strelka"
                    ], vcf]
                }

    strelka_vcf_tbi = Channel.empty().mix(
                        MERGE_STRELKA.out.tbi,
                        strelka_vcf_tbi.no_intervals)
                    .map{ meta, tbi ->
                        [[
                            id:             meta.sample,
                            num_intervals:  meta.num_intervals,
                            patient:        meta.patient,
                            sample:         meta.sample,
                            sex:            meta.sex,
                            status:         meta.status,
                            variantcaller:  "strelka"
                        ], tbi]
                    }

    ch_versions = ch_versions.mix(MERGE_STRELKA.out.versions)
    ch_versions = ch_versions.mix(MERGE_STRELKA_GENOME.out.versions)
    ch_versions = ch_versions.mix(STRELKA_SINGLE.out.versions)

    emit:
    strelka_vcf
    strelka_vcf_tbi
    versions = ch_versions
}
