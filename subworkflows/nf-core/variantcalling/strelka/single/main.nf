include { GATK4_MERGEVCFS as MERGE_STRELKA        } from '../../../../../modules/nf-core/modules/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_STRELKA_GENOME } from '../../../../../modules/nf-core/modules/gatk4/mergevcfs/main'
include { STRELKA_GERMLINE as STRELKA_SINGLE      } from '../../../../../modules/nf-core/modules/strelka/germline/main'

workflow RUN_STRELKA_SINGLE {
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

    STRELKA_SINGLE.out.genome_vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{strelka_genome_vcf}

    MERGE_STRELKA(
        strelka_vcf.intervals
            .map{ meta, vcf ->
                new_meta = [patient:meta.patient, sample:meta.sample, status:meta.status, gender:meta.gender, id:meta.sample, num_intervals:meta.num_intervals]

                [groupKey(new_meta, meta.num_intervals), vcf]
            }.groupTuple(),
        dict
    )

    MERGE_STRELKA_GENOME(
        strelka_genome_vcf.intervals
            .map{ meta, vcf ->

                [groupKey([patient:meta.patient, sample:meta.sample, status:meta.status, gender:meta.gender, id:meta.sample, num_intervals:meta.num_intervals],
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
                    [[patient:meta.patient, sample:meta.sample, status:meta.status, gender:meta.gender, id:meta.sample, num_intervals:meta.num_intervals, variantcaller:"Strelka"], vcf]
                }

    ch_versions = ch_versions.mix(MERGE_STRELKA.out.versions)
    ch_versions = ch_versions.mix(MERGE_STRELKA_GENOME.out.versions)
    ch_versions = ch_versions.mix(STRELKA_SINGLE.out.versions)

    emit:
    strelka_vcf
    versions = ch_versions
}
