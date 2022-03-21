include { BGZIP as BGZIP_VC_STRELKA                 } from '../../modules/local/bgzip'
include { BGZIP as BGZIP_VC_STRELKA_GENOME          } from '../../modules/local/bgzip'
include { CONCAT_VCF as CONCAT_STRELKA              } from '../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_STRELKA_GENOME       } from '../../modules/local/concat_vcf/main'
include { STRELKA_GERMLINE                          } from '../../modules/nf-core/modules/strelka/germline/main'

workflow RUN_STRELKA {
    take:
    cram_recalibrated_intervals_gz_tbi
    fasta
    fasta_fai
    num_intervals

    main:

    ch_versions = Channel.empty()

    // TODO: Research if splitting by intervals is ok, we pretend for now it is fine.
    // Seems to be the consensus on upstream modules implementation too

    STRELKA_GERMLINE(
        cram_recalibrated_intervals_gz_tbi,
        fasta,
        fasta_fai)

    // Figure out if using intervals or no_intervals
    STRELKA_GERMLINE.out.vcf.groupTuple(size: num_intervals)
        .branch{
            intervals:    it[1].size() > 1
            no_intervals: it[1].size() == 1
        }.set{strelka_vcf}

    STRELKA_GERMLINE.out.genome_vcf.groupTuple(size: num_intervals)
        .branch{
            intervals:    it[1].size() > 1
            no_intervals: it[1].size() == 1
        }.set{strelka_genome_vcf}

    // Only when using intervals
    BGZIP_VC_STRELKA(STRELKA_GERMLINE.out.vcf)

    CONCAT_STRELKA(
        BGZIP_VC_STRELKA.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_combine_gz)

    BGZIP_VC_STRELKA_GENOME(STRELKA_GERMLINE.out.genome_vcf)

    CONCAT_STRELKA_GENOME(
        BGZIP_VC_STRELKA_GENOME.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_combine_gz)

    strelka_vcf = Channel.empty().mix(
        CONCAT_STRELKA.out.vcf,
        CONCAT_STRELKA_GENOME.out.vcf,
        strelka_genome_vcf.no_intervals,

    ch_versions = ch_versions.mix(BGZIP_VC_STRELKA.out.versions)
    ch_versions = ch_versions.mix(CONCAT_STRELKA.out.versions)
    ch_versions = ch_versions.mix(STRELKA_GERMLINE.out.versions)

    emit:
    versions = ch_versions
    strelka_vcf
}
