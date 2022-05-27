include { TABIX_BGZIP as BGZIP_VC_STRELKA        } from '../../../../../modules/nf-core/modules/tabix/bgzip/main'
include { TABIX_BGZIP as BGZIP_VC_STRELKA_GENOME } from '../../../../../modules/nf-core/modules/tabix/bgzip/main'
include { CONCAT_VCF as CONCAT_STRELKA           } from '../../../../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_STRELKA_GENOME    } from '../../../../../modules/local/concat_vcf/main'
include { STRELKA_GERMLINE as STRELKA_SINGLE     } from '../../../../../modules/nf-core/modules/strelka/germline/main'

workflow RUN_STRELKA_SINGLE {
    take:
    cram                     // channel: [mandatory] [meta, cram, crai, interval.bed.gz, interval.bed.gz.tbi]
    fasta                    // channel: [mandatory]
    fasta_fai                // channel: [mandatory]
    intervals_bed_gz         // channel: [optional]  Contains a bed.gz file of all intervals combined provided with the cram input(s). Mandatory if interval files are used.

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

    // Only when using intervals
    BGZIP_VC_STRELKA(strelka_vcf.intervals)

    CONCAT_STRELKA(
        BGZIP_VC_STRELKA.out.output
            .map{ meta, vcf ->

                [groupKey([patient:meta.patient, sample:meta.sample, status:meta.status, gender:meta.gender, id:meta.sample, num_intervals:meta.num_intervals],
                        meta.num_intervals),
                vcf]

            }.groupTuple(),
        fasta_fai,
        intervals_bed_gz)

    BGZIP_VC_STRELKA_GENOME(strelka_genome_vcf.intervals)

    CONCAT_STRELKA_GENOME(
        BGZIP_VC_STRELKA_GENOME.out.output
            .map{ meta, vcf ->
                [groupKey([patient:meta.patient, sample:meta.sample, status:meta.status, gender:meta.gender, id:meta.sample, num_intervals:meta.num_intervals],
                        meta.num_intervals),
                vcf]

            }.groupTuple(),
        fasta_fai,
        intervals_bed_gz)

    // Mix output channels for "no intervals" and "with intervals" results
    // Only strelka variant vcf should get annotated
    strelka_vcf = Channel.empty().mix(
                    CONCAT_STRELKA.out.vcf,
                    strelka_vcf.no_intervals)
                .map{ meta, vcf ->
                    [[patient:meta.patient, sample:meta.sample, status:meta.status, gender:meta.gender, id:meta.sample, num_intervals:meta.num_intervals, variantcaller:"Strelka"], vcf]
                }

    ch_versions = ch_versions.mix(BGZIP_VC_STRELKA.out.versions)
    ch_versions = ch_versions.mix(BGZIP_VC_STRELKA_GENOME.out.versions)
    ch_versions = ch_versions.mix(CONCAT_STRELKA.out.versions)
    ch_versions = ch_versions.mix(CONCAT_STRELKA_GENOME.out.versions)
    ch_versions = ch_versions.mix(STRELKA_SINGLE.out.versions)

    emit:
    strelka_vcf
    versions = ch_versions
}
