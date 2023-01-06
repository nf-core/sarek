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
    versions = Channel.empty()

    STRELKA_SINGLE(cram, fasta, fasta_fai)

    // Figure out if using intervals or no_intervals
    vcf = STRELKA_SINGLE.out.vcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    genome_vcf = STRELKA_SINGLE.out.genome_vcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    MERGE_STRELKA(
        vcf.intervals
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
        dict.map{ it -> [ [ id:'dict' ], it ] })

    MERGE_STRELKA_GENOME(
        genome_vcf.intervals
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
        dict.map{ it -> [ [ id:'dict' ], it ] })

    // Mix output channels for "no intervals" and "with intervals" results
    // Only strelka variant vcf should get annotated
    vcf = Channel.empty().mix(
                    MERGE_STRELKA.out.vcf,
                    vcf.no_intervals)
                .map{ meta, vcf ->
                    [ [
                        id:             meta.sample,
                        num_intervals:  meta.num_intervals,
                        patient:        meta.patient,
                        sample:         meta.sample,
                        sex:            meta.sex,
                        status:         meta.status,
                        variantcaller:  "strelka"
                    ], vcf ]
                }

    versions = versions.mix(MERGE_STRELKA.out.versions)
    versions = versions.mix(MERGE_STRELKA_GENOME.out.versions)
    versions = versions.mix(STRELKA_SINGLE.out.versions)

    emit:
    vcf

    versions
}
