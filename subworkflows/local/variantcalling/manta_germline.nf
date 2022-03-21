include { BGZIP as BGZIP_VC_MANTA_DIPLOID           } from '../../../modules/local/bgzip'
include { BGZIP as BGZIP_VC_MANTA_SMALL_INDELS      } from '../../../modules/local/bgzip'
include { BGZIP as BGZIP_VC_MANTA_SV                } from '../../../modules/local/bgzip'
include { CONCAT_VCF as CONCAT_MANTA_DIPLOID        } from '../../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_MANTA_SMALL_INDELS   } from '../../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_MANTA_SV             } from '../../../modules/local/concat_vcf/main'
include { MANTA_GERMLINE                            } from '../../../modules/local/manta/germline/main'

workflow RUN_MANTA {
    take:
    cram_recalibrated_intervals_gz_tbi
    fasta
    fasta_fai
    num_intervals
    intervals_bed_combine_gz

    main:

    ch_versions = Channel.empty()
    // TODO: Research if splitting by intervals is ok, we pretend for now it is fine.
    // Seems to be the consensus on upstream modules implementation too

    MANTA_GERMLINE(
        cram_recalibrated_intervals_gz_tbi,
        fasta,
        fasta_fai)

    // Figure out if using intervals or no_intervals
    MANTA_GERMLINE.out.candidate_small_indels_vcf.groupTuple(size: num_intervals)
        .branch{
            intervals:    it[1].size() > 1
            no_intervals: it[1].size() == 1
        }.set{manta_small_indels_vcf}

    MANTA_GERMLINE.out.candidate_sv_vcf.groupTuple(size: num_intervals)
        .branch{
            intervals:    it[1].size() > 1
            no_intervals: it[1].size() == 1
        }.set{manta_sv_vcf}

    MANTA_GERMLINE.out.diploid_sv_vcf.groupTuple(size: num_intervals)
        .branch{
            intervals:    it[1].size() > 1
            no_intervals: it[1].size() == 1
        }.set{manta_diploid_sv_vcf}

    // Only when using intervals
    BGZIP_VC_MANTA_DIPLOID(MANTA_GERMLINE.out.diploid_sv_vcf)

    CONCAT_MANTA_DIPLOID(
        BGZIP_VC_MANTA_DIPLOID.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_combine_gz)

    BGZIP_VC_MANTA_SMALL_INDELS(MANTA_GERMLINE.out.candidate_small_indels_vcf)

    CONCAT_MANTA_SMALL_INDELS(
        BGZIP_VC_MANTA_SMALL_INDELS.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_combine_gz)

    BGZIP_VC_MANTA_SV(MANTA_GERMLINE.out.candidate_sv_vcf)

    CONCAT_MANTA_SV(
        BGZIP_VC_MANTA_SV.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_combine_gz)

    manta_vcf = Channel.empty().mix(
        CONCAT_MANTA_DIPLOID.out.vcf,
        CONCAT_MANTA_SMALL_INDELS.out.vcf,
        CONCAT_MANTA_SV.out.vcf,
        manta_diploid_sv_vcf.no_intervals,
        manta_small_indels_vcf.no_intervals,
        manta_sv_vcf.no_intervals)

    ch_versions = ch_versions.mix(BGZIP_VC_MANTA_DIPLOID.out.versions)
    ch_versions = ch_versions.mix(BGZIP_VC_MANTA_SMALL_INDELS.out.versions)
    ch_versions = ch_versions.mix(BGZIP_VC_MANTA_SV.out.versions)
    ch_versions = ch_versions.mix(CONCAT_MANTA_DIPLOID.out.versions)
    ch_versions = ch_versions.mix(CONCAT_MANTA_SMALL_INDELS.out.versions)
    ch_versions = ch_versions.mix(CONCAT_MANTA_SV.out.versions)
    ch_versions = ch_versions.mix(MANTA_GERMLINE.out.versions)

    emit:
    versions = ch_versions
    manta_vcf
}
