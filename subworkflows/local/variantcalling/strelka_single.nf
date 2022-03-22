include { BGZIP as BGZIP_VC_STRELKA                 } from '../../../modules/local/bgzip'
include { BGZIP as BGZIP_VC_STRELKA_GENOME          } from '../../../modules/local/bgzip'
include { CONCAT_VCF as CONCAT_STRELKA        } from '../../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_STRELKA_GENOME } from '../../../modules/local/concat_vcf/main'
include { STRELKA_GERMLINE                    } from '../../../modules/nf-core/modules/strelka/germline/main'

// TODO: Research if splitting by intervals is ok, we pretend for now it is fine.
// Seems to be the consensus on upstream modules implementation too
workflow RUN_STRELKA_SINGLE {
    take:
    cram                     // channel: [mandatory] [meta, cram, crai, interval.bed.gz, interval.bed.gz.tbi]
    fasta                    // channel: [mandatory]
    fasta_fai                // channel: [mandatory]
    intervals_bed_gz         // channel: [optional]  Contains a bed.gz file of all intervals combined provided with the cram input(s). Mandatory if interval files are used.
    num_intervals            //     val: [optional]  Number of used intervals, mandatory when intervals are provided.

    main:

    ch_versions = Channel.empty()

    STRELKA_GERMLINE(cram, fasta, fasta_fai)

    // Figure out if using intervals or no_intervals
    STRELKA_GERMLINE.out.vcf.branch{
            intervals:    num_intervals > 1
            no_intervals: num_intervals == 1
        }
        .set{strelka_vcf}

    STRELKA_GERMLINE.out.genome_vcf.branch{
            intervals:    num_intervals > 1
            no_intervals: num_intervals == 1
        }.set{strelka_genome_vcf}

    // Only when using intervals
    BGZIP_VC_STRELKA(strelka_vcf.intervals)

    CONCAT_STRELKA(
        BGZIP_VC_STRELKA.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_gz)

    BGZIP_VC_STRELKA_GENOME(strelka_genome_vcf.intervals)

    CONCAT_STRELKA_GENOME(
        BGZIP_VC_STRELKA_GENOME.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_gz)

    // Mix output channels for "no intervals" and "with intervals" results
    strelka_vcf = Channel.empty().mix(
        CONCAT_STRELKA.out.vcf,
        CONCAT_STRELKA_GENOME.out.vcf,
        strelka_genome_vcf.no_intervals,
        strelka_vcf.no_intervals)

    ch_versions = ch_versions.mix(BGZIP_VC_STRELKA.out.versions)
    ch_versions = ch_versions.mix(CONCAT_STRELKA.out.versions)
    ch_versions = ch_versions.mix(STRELKA_GERMLINE.out.versions)

    emit:
    strelka_vcf
    versions = ch_versions
}
